% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPPRR4 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:14
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RPPRR4_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPRR4_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:14:30
	% EndTime: 2020-11-04 20:14:30
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:14:30
	% EndTime: 2020-11-04 20:14:30
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t43 = cos(qJ(1));
	t42 = sin(qJ(1));
	t1 = [0, 0, 1, pkin(5) + 0; t42, t43, 0, 0; -t43, t42, 0, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:14:30
	% EndTime: 2020-11-04 20:14:30
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->9), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t47 = cos(qJ(1));
	t46 = sin(qJ(1));
	t45 = cos(pkin(8));
	t44 = sin(pkin(8));
	t1 = [t44, t45, 0, pkin(5) + 0; t46 * t45, -t46 * t44, -t47, t46 * pkin(1) - t47 * qJ(2) + 0; -t47 * t45, t47 * t44, -t46, -t47 * pkin(1) - t46 * qJ(2) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:14:30
	% EndTime: 2020-11-04 20:14:30
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (18->16), mult. (26->20), div. (0->0), fcn. (39->6), ass. (0->12)
	t49 = sin(pkin(9));
	t53 = sin(qJ(1));
	t58 = t53 * t49;
	t51 = cos(pkin(9));
	t57 = t53 * t51;
	t54 = cos(qJ(1));
	t56 = t54 * t49;
	t55 = t54 * t51;
	t52 = cos(pkin(8));
	t50 = sin(pkin(8));
	t48 = pkin(2) * t52 + t50 * qJ(3) + pkin(1);
	t1 = [t50 * t51, -t50 * t49, -t52, t50 * pkin(2) - t52 * qJ(3) + pkin(5) + 0; t52 * t57 - t56, -t52 * t58 - t55, t53 * t50, -t54 * qJ(2) + t48 * t53 + 0; -t52 * t55 - t58, t52 * t56 - t57, -t54 * t50, -t53 * qJ(2) - t48 * t54 + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:14:30
	% EndTime: 2020-11-04 20:14:30
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (36->20), mult. (31->22), div. (0->0), fcn. (44->8), ass. (0->16)
	t63 = pkin(9) + qJ(4);
	t61 = sin(t63);
	t67 = sin(qJ(1));
	t73 = t67 * t61;
	t62 = cos(t63);
	t72 = t67 * t62;
	t68 = cos(qJ(1));
	t71 = t68 * t61;
	t70 = t68 * t62;
	t60 = cos(pkin(9)) * pkin(3) + pkin(2);
	t64 = sin(pkin(8));
	t65 = cos(pkin(8));
	t66 = qJ(3) + pkin(6);
	t69 = t60 * t65 + t66 * t64 + pkin(1);
	t59 = -sin(pkin(9)) * pkin(3) - qJ(2);
	t1 = [t64 * t62, -t64 * t61, -t65, t64 * t60 - t65 * t66 + pkin(5) + 0; t65 * t72 - t71, -t65 * t73 - t70, t67 * t64, t59 * t68 + t69 * t67 + 0; -t65 * t70 - t73, t65 * t71 - t72, -t68 * t64, t59 * t67 - t69 * t68 + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:14:30
	% EndTime: 2020-11-04 20:14:30
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (59->23), mult. (42->24), div. (0->0), fcn. (55->10), ass. (0->17)
	t80 = pkin(9) + qJ(4);
	t78 = qJ(5) + t80;
	t76 = sin(t78);
	t83 = sin(qJ(1));
	t90 = t83 * t76;
	t77 = cos(t78);
	t89 = t83 * t77;
	t84 = cos(qJ(1));
	t88 = t84 * t76;
	t87 = t84 * t77;
	t86 = -qJ(2) - pkin(4) * sin(t80) - sin(pkin(9)) * pkin(3);
	t74 = pkin(4) * cos(t80) + cos(pkin(9)) * pkin(3) + pkin(2);
	t79 = -pkin(7) - pkin(6) - qJ(3);
	t81 = sin(pkin(8));
	t82 = cos(pkin(8));
	t85 = t74 * t82 - t79 * t81 + pkin(1);
	t1 = [t81 * t77, -t81 * t76, -t82, t81 * t74 + t82 * t79 + pkin(5) + 0; t82 * t89 - t88, -t82 * t90 - t87, t83 * t81, t85 * t83 + t86 * t84 + 0; -t82 * t87 - t90, t82 * t88 - t89, -t84 * t81, t86 * t83 - t85 * t84 + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end