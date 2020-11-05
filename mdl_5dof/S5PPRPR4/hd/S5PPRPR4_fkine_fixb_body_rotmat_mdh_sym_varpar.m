% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5PPRPR4 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 19:53
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5PPRPR4_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR4_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PPRPR4_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRPR4_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:53:22
	% EndTime: 2020-11-04 19:53:22
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:53:22
	% EndTime: 2020-11-04 19:53:22
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t38 = cos(pkin(7));
	t37 = sin(pkin(7));
	t1 = [t38, -t37, 0, 0; t37, t38, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:53:22
	% EndTime: 2020-11-04 19:53:22
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (6->6), mult. (4->4), div. (0->0), fcn. (8->2), ass. (0->3)
	t40 = cos(pkin(7));
	t39 = sin(pkin(7));
	t1 = [t40, 0, t39, t40 * pkin(1) + t39 * qJ(2) + 0; t39, 0, -t40, t39 * pkin(1) - t40 * qJ(2) + 0; 0, 1, 0, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:53:22
	% EndTime: 2020-11-04 19:53:22
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (12->10), mult. (12->8), div. (0->0), fcn. (20->4), ass. (0->8)
	t47 = pkin(1) + pkin(2);
	t46 = cos(qJ(3));
	t45 = sin(qJ(3));
	t44 = cos(pkin(7));
	t43 = sin(pkin(7));
	t42 = t43 * t46 - t44 * t45;
	t41 = -t43 * t45 - t44 * t46;
	t1 = [-t41, t42, 0, t43 * qJ(2) + t47 * t44 + 0; t42, t41, 0, -t44 * qJ(2) + t47 * t43 + 0; 0, 0, -1, -pkin(5) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:53:22
	% EndTime: 2020-11-04 19:53:22
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (25->19), mult. (32->20), div. (0->0), fcn. (46->6), ass. (0->12)
	t58 = pkin(1) + pkin(2);
	t57 = cos(qJ(3));
	t56 = sin(qJ(3));
	t55 = cos(pkin(7));
	t54 = cos(pkin(8));
	t53 = sin(pkin(7));
	t52 = sin(pkin(8));
	t51 = t55 * pkin(3) - qJ(4) * t53;
	t50 = pkin(3) * t53 + t55 * qJ(4);
	t49 = t53 * t56 + t55 * t57;
	t48 = -t53 * t57 + t55 * t56;
	t1 = [t49 * t54, -t49 * t52, t48, t53 * qJ(2) + t50 * t56 + t51 * t57 + t58 * t55 + 0; -t48 * t54, t48 * t52, t49, -t55 * qJ(2) + t50 * t57 - t51 * t56 + t58 * t53 + 0; -t52, -t54, 0, -pkin(5) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:53:22
	% EndTime: 2020-11-04 19:53:22
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (37->21), mult. (35->18), div. (0->0), fcn. (53->8), ass. (0->13)
	t71 = cos(qJ(3));
	t70 = sin(pkin(7));
	t69 = pkin(1) + pkin(2);
	t68 = sin(qJ(3));
	t67 = -pkin(6) - qJ(4);
	t66 = cos(pkin(7));
	t65 = pkin(8) + qJ(5);
	t64 = cos(t65);
	t63 = sin(t65);
	t62 = cos(pkin(8)) * pkin(4) + pkin(3);
	t60 = t66 * t71 + t70 * t68;
	t59 = t66 * t68 - t70 * t71;
	t1 = [t60 * t64, -t60 * t63, t59, t70 * qJ(2) - t59 * t67 + t60 * t62 + t69 * t66 + 0; -t59 * t64, t59 * t63, t60, -t66 * qJ(2) - t59 * t62 - t60 * t67 + t69 * t70 + 0; -t63, -t64, 0, -sin(pkin(8)) * pkin(4) - pkin(5) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end