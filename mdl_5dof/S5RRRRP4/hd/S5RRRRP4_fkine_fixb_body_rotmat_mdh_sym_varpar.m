% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRRRP4 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:45
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RRRRP4_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP4_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRP4_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP4_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:45:51
	% EndTime: 2020-11-04 20:45:51
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:45:51
	% EndTime: 2020-11-04 20:45:51
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t40 = cos(qJ(1));
	t39 = sin(qJ(1));
	t1 = [t40, -t39, 0, 0; t39, t40, 0, 0; 0, 0, 1, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:45:51
	% EndTime: 2020-11-04 20:45:51
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->6), mult. (2->2), div. (0->0), fcn. (6->4), ass. (0->4)
	t43 = qJ(1) + qJ(2);
	t42 = cos(t43);
	t41 = sin(t43);
	t1 = [t42, -t41, 0, cos(qJ(1)) * pkin(1) + 0; t41, t42, 0, sin(qJ(1)) * pkin(1) + 0; 0, 0, 1, pkin(6) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:45:51
	% EndTime: 2020-11-04 20:45:51
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (21->12), mult. (10->10), div. (0->0), fcn. (18->6), ass. (0->6)
	t48 = cos(qJ(3));
	t47 = sin(qJ(3));
	t46 = qJ(1) + qJ(2);
	t45 = cos(t46);
	t44 = sin(t46);
	t1 = [t45 * t48, -t45 * t47, t44, t45 * pkin(2) + t44 * pkin(7) + cos(qJ(1)) * pkin(1) + 0; t44 * t48, -t44 * t47, -t45, t44 * pkin(2) - t45 * pkin(7) + sin(qJ(1)) * pkin(1) + 0; t47, t48, 0, pkin(6) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:45:51
	% EndTime: 2020-11-04 20:45:51
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (32->16), mult. (13->12), div. (0->0), fcn. (21->8), ass. (0->9)
	t56 = -pkin(8) - pkin(7);
	t55 = qJ(1) + qJ(2);
	t54 = qJ(3) + qJ(4);
	t53 = cos(t55);
	t52 = cos(t54);
	t51 = sin(t55);
	t50 = sin(t54);
	t49 = cos(qJ(3)) * pkin(3) + pkin(2);
	t1 = [t53 * t52, -t53 * t50, t51, t53 * t49 - t51 * t56 + cos(qJ(1)) * pkin(1) + 0; t51 * t52, -t51 * t50, -t53, t51 * t49 + t53 * t56 + sin(qJ(1)) * pkin(1) + 0; t50, t52, 0, sin(qJ(3)) * pkin(3) + pkin(6) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:45:51
	% EndTime: 2020-11-04 20:45:51
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (47->19), mult. (23->16), div. (0->0), fcn. (31->8), ass. (0->9)
	t62 = qJ(3) + qJ(4);
	t58 = sin(t62);
	t60 = cos(t62);
	t65 = pkin(4) * t60 + qJ(5) * t58 + cos(qJ(3)) * pkin(3) + pkin(2);
	t64 = -pkin(8) - pkin(7);
	t63 = qJ(1) + qJ(2);
	t61 = cos(t63);
	t59 = sin(t63);
	t1 = [t61 * t60, t59, t61 * t58, cos(qJ(1)) * pkin(1) - t59 * t64 + 0 + t65 * t61; t59 * t60, -t61, t59 * t58, sin(qJ(1)) * pkin(1) + t61 * t64 + 0 + t65 * t59; t58, 0, -t60, t58 * pkin(4) - t60 * qJ(5) + sin(qJ(3)) * pkin(3) + pkin(6) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end