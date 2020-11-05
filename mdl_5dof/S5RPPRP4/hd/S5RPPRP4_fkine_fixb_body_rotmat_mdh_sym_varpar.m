% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPPRP4 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:12
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RPPRP4_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP4_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPRP4_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP4_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [7x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:12:53
	% EndTime: 2020-11-04 20:12:53
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:12:53
	% EndTime: 2020-11-04 20:12:53
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t39 = cos(qJ(1));
	t38 = sin(qJ(1));
	t1 = [t39, -t38, 0, 0; t38, t39, 0, 0; 0, 0, 1, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:12:53
	% EndTime: 2020-11-04 20:12:53
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (6->6), mult. (4->4), div. (0->0), fcn. (8->2), ass. (0->3)
	t41 = cos(qJ(1));
	t40 = sin(qJ(1));
	t1 = [t41, 0, t40, t41 * pkin(1) + t40 * qJ(2) + 0; t40, 0, -t41, t40 * pkin(1) - t41 * qJ(2) + 0; 0, 1, 0, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:12:53
	% EndTime: 2020-11-04 20:12:53
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (12->10), mult. (12->8), div. (0->0), fcn. (20->4), ass. (0->8)
	t48 = pkin(1) + pkin(2);
	t47 = cos(qJ(1));
	t46 = sin(qJ(1));
	t45 = cos(pkin(7));
	t44 = sin(pkin(7));
	t43 = -t47 * t44 + t46 * t45;
	t42 = -t46 * t44 - t47 * t45;
	t1 = [-t42, t43, 0, t46 * qJ(2) + t48 * t47 + 0; t43, t42, 0, -t47 * qJ(2) + t48 * t46 + 0; 0, 0, -1, -qJ(3) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:12:53
	% EndTime: 2020-11-04 20:12:54
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (26->17), mult. (28->16), div. (0->0), fcn. (42->6), ass. (0->11)
	t59 = cos(qJ(1));
	t58 = cos(qJ(4));
	t57 = sin(qJ(1));
	t56 = sin(qJ(4));
	t55 = cos(pkin(7));
	t54 = sin(pkin(7));
	t52 = -t54 * pkin(3) + t55 * pkin(6) - qJ(2);
	t51 = t55 * pkin(3) + t54 * pkin(6) + pkin(1) + pkin(2);
	t50 = t57 * t54 + t59 * t55;
	t49 = t59 * t54 - t57 * t55;
	t1 = [t50 * t58, -t50 * t56, t49, t51 * t59 - t52 * t57 + 0; -t49 * t58, t49 * t56, t50, t51 * t57 + t52 * t59 + 0; -t56, -t58, 0, -qJ(3) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:12:54
	% EndTime: 2020-11-04 20:12:54
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (35->20), mult. (37->18), div. (0->0), fcn. (51->6), ass. (0->13)
	t64 = sin(pkin(7));
	t65 = cos(pkin(7));
	t66 = qJ(5) + pkin(6);
	t69 = cos(qJ(4));
	t72 = pkin(4) * t69 + pkin(3);
	t73 = t72 * t64 - t66 * t65 + qJ(2);
	t70 = cos(qJ(1));
	t68 = sin(qJ(1));
	t67 = sin(qJ(4));
	t62 = t68 * t64 + t70 * t65;
	t61 = t70 * t64 - t68 * t65;
	t60 = t66 * t64 + t72 * t65 + pkin(1) + pkin(2);
	t1 = [t62 * t69, -t62 * t67, t61, t60 * t70 + t73 * t68 + 0; -t61 * t69, t61 * t67, t62, t60 * t68 - t73 * t70 + 0; -t67, -t69, 0, -t67 * pkin(4) + pkin(5) - qJ(3) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end