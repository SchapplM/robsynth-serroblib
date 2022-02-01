% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPRRP3 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for the body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:30
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RPRRP3_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP3_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRRP3_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP3_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:29:27
	% EndTime: 2022-01-23 09:29:27
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:29:27
	% EndTime: 2022-01-23 09:29:27
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t37 = cos(qJ(1));
	t36 = sin(qJ(1));
	t1 = [t37, -t36, 0, 0; t36, t37, 0, 0; 0, 0, 1, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:29:27
	% EndTime: 2022-01-23 09:29:27
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (9->6), mult. (2->2), div. (0->0), fcn. (6->4), ass. (0->4)
	t40 = qJ(1) + pkin(8);
	t39 = cos(t40);
	t38 = sin(t40);
	t1 = [t39, -t38, 0, cos(qJ(1)) * pkin(1) + 0; t38, t39, 0, sin(qJ(1)) * pkin(1) + 0; 0, 0, 1, qJ(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:29:27
	% EndTime: 2022-01-23 09:29:27
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (21->12), mult. (10->10), div. (0->0), fcn. (18->6), ass. (0->6)
	t45 = cos(qJ(3));
	t44 = sin(qJ(3));
	t43 = qJ(1) + pkin(8);
	t42 = cos(t43);
	t41 = sin(t43);
	t1 = [t42 * t45, -t42 * t44, t41, t42 * pkin(2) + t41 * pkin(6) + cos(qJ(1)) * pkin(1) + 0; t41 * t45, -t41 * t44, -t42, t41 * pkin(2) - t42 * pkin(6) + sin(qJ(1)) * pkin(1) + 0; t44, t45, 0, qJ(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:29:27
	% EndTime: 2022-01-23 09:29:27
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (32->16), mult. (13->12), div. (0->0), fcn. (21->8), ass. (0->9)
	t53 = -pkin(7) - pkin(6);
	t52 = qJ(3) + qJ(4);
	t51 = qJ(1) + pkin(8);
	t50 = cos(t52);
	t49 = sin(t52);
	t48 = cos(t51);
	t47 = sin(t51);
	t46 = cos(qJ(3)) * pkin(3) + pkin(2);
	t1 = [t48 * t50, -t48 * t49, t47, t48 * t46 - t47 * t53 + cos(qJ(1)) * pkin(1) + 0; t47 * t50, -t47 * t49, -t48, t47 * t46 + t48 * t53 + sin(qJ(1)) * pkin(1) + 0; t49, t50, 0, sin(qJ(3)) * pkin(3) + qJ(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:29:27
	% EndTime: 2022-01-23 09:29:27
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (40->18), mult. (16->14), div. (0->0), fcn. (24->8), ass. (0->9)
	t61 = qJ(3) + qJ(4);
	t60 = qJ(1) + pkin(8);
	t59 = -qJ(5) - pkin(7) - pkin(6);
	t58 = cos(t61);
	t57 = sin(t61);
	t56 = cos(t60);
	t55 = sin(t60);
	t54 = pkin(4) * t58 + cos(qJ(3)) * pkin(3) + pkin(2);
	t1 = [t56 * t58, -t56 * t57, t55, t56 * t54 - t55 * t59 + cos(qJ(1)) * pkin(1) + 0; t55 * t58, -t55 * t57, -t56, t55 * t54 + t56 * t59 + sin(qJ(1)) * pkin(1) + 0; t57, t58, 0, pkin(4) * t57 + sin(qJ(3)) * pkin(3) + qJ(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end