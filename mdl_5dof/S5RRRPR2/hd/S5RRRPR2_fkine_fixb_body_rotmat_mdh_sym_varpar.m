% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRRPR2 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for the body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 11:31
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RRRPR2_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR2_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRPR2_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR2_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 11:30:27
	% EndTime: 2022-01-20 11:30:27
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 11:30:27
	% EndTime: 2022-01-20 11:30:27
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t41 = cos(qJ(1));
	t40 = sin(qJ(1));
	t1 = [t41, -t40, 0, 0; t40, t41, 0, 0; 0, 0, 1, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 11:30:27
	% EndTime: 2022-01-20 11:30:27
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (9->6), mult. (2->2), div. (0->0), fcn. (6->4), ass. (0->4)
	t44 = qJ(1) + qJ(2);
	t43 = cos(t44);
	t42 = sin(t44);
	t1 = [t43, -t42, 0, cos(qJ(1)) * pkin(1) + 0; t42, t43, 0, sin(qJ(1)) * pkin(1) + 0; 0, 0, 1, pkin(6) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 11:30:27
	% EndTime: 2022-01-20 11:30:27
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (18->10), mult. (4->4), div. (0->0), fcn. (8->6), ass. (0->5)
	t48 = qJ(1) + qJ(2);
	t47 = qJ(3) + t48;
	t46 = cos(t47);
	t45 = sin(t47);
	t1 = [t46, -t45, 0, pkin(2) * cos(t48) + cos(qJ(1)) * pkin(1) + 0; t45, t46, 0, pkin(2) * sin(t48) + sin(qJ(1)) * pkin(1) + 0; 0, 0, 1, pkin(7) + pkin(6) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 11:30:27
	% EndTime: 2022-01-20 11:30:27
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (29->14), mult. (6->6), div. (0->0), fcn. (10->8), ass. (0->6)
	t53 = qJ(1) + qJ(2);
	t52 = qJ(3) + t53;
	t51 = pkin(9) + t52;
	t50 = cos(t51);
	t49 = sin(t51);
	t1 = [t50, -t49, 0, pkin(3) * cos(t52) + pkin(2) * cos(t53) + cos(qJ(1)) * pkin(1) + 0; t49, t50, 0, pkin(3) * sin(t52) + pkin(2) * sin(t53) + sin(qJ(1)) * pkin(1) + 0; 0, 0, 1, qJ(4) + pkin(7) + pkin(6) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 11:30:27
	% EndTime: 2022-01-20 11:30:27
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (53->20), mult. (14->14), div. (0->0), fcn. (22->10), ass. (0->8)
	t58 = qJ(1) + qJ(2);
	t57 = qJ(3) + t58;
	t60 = cos(qJ(5));
	t59 = sin(qJ(5));
	t56 = pkin(9) + t57;
	t55 = cos(t56);
	t54 = sin(t56);
	t1 = [t55 * t60, -t55 * t59, t54, t55 * pkin(4) + t54 * pkin(8) + pkin(3) * cos(t57) + pkin(2) * cos(t58) + cos(qJ(1)) * pkin(1) + 0; t54 * t60, -t54 * t59, -t55, t54 * pkin(4) - t55 * pkin(8) + pkin(3) * sin(t57) + pkin(2) * sin(t58) + sin(qJ(1)) * pkin(1) + 0; t59, t60, 0, qJ(4) + pkin(7) + pkin(6) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end