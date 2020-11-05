% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPPPRR4 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:23
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RPPPRR4_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR4_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPPRR4_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR4_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:23:58
	% EndTime: 2020-11-04 21:23:58
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:23:58
	% EndTime: 2020-11-04 21:23:58
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t53 = cos(qJ(1));
	t52 = sin(qJ(1));
	t1 = [t53, -t52, 0, 0; t52, t53, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:23:58
	% EndTime: 2020-11-04 21:23:58
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (6->6), mult. (4->4), div. (0->0), fcn. (8->2), ass. (0->3)
	t55 = cos(qJ(1));
	t54 = sin(qJ(1));
	t1 = [t55, 0, t54, pkin(1) * t55 + qJ(2) * t54 + 0; t54, 0, -t55, pkin(1) * t54 - qJ(2) * t55 + 0; 0, 1, 0, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:23:58
	% EndTime: 2020-11-04 21:23:58
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (12->10), mult. (12->8), div. (0->0), fcn. (20->4), ass. (0->8)
	t62 = pkin(1) + pkin(2);
	t61 = cos(qJ(1));
	t60 = sin(qJ(1));
	t59 = cos(pkin(9));
	t58 = sin(pkin(9));
	t57 = -t61 * t58 + t60 * t59;
	t56 = -t60 * t58 - t61 * t59;
	t1 = [-t56, t57, 0, t60 * qJ(2) + t62 * t61 + 0; t57, t56, 0, -t61 * qJ(2) + t62 * t60 + 0; 0, 0, -1, -qJ(3) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:23:58
	% EndTime: 2020-11-04 21:23:58
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (20->14), mult. (20->12), div. (0->0), fcn. (28->4), ass. (0->9)
	t70 = cos(qJ(1));
	t69 = sin(qJ(1));
	t68 = cos(pkin(9));
	t67 = sin(pkin(9));
	t66 = -t67 * pkin(3) + t68 * qJ(4) - qJ(2);
	t65 = t68 * pkin(3) + t67 * qJ(4) + pkin(1) + pkin(2);
	t64 = t70 * t67 - t69 * t68;
	t63 = -t69 * t67 - t70 * t68;
	t1 = [0, t63, t64, t65 * t70 - t66 * t69 + 0; 0, t64, -t63, t65 * t69 + t66 * t70 + 0; -1, 0, 0, -qJ(3) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:23:58
	% EndTime: 2020-11-04 21:23:58
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (28->17), mult. (28->16), div. (0->0), fcn. (42->6), ass. (0->12)
	t74 = sin(pkin(9));
	t75 = cos(pkin(9));
	t80 = pkin(3) + pkin(7);
	t81 = t75 * qJ(4) - t80 * t74 - qJ(2);
	t79 = cos(qJ(1));
	t78 = cos(qJ(5));
	t77 = sin(qJ(1));
	t76 = sin(qJ(5));
	t73 = t77 * t74 + t79 * t75;
	t72 = t79 * t74 - t77 * t75;
	t71 = t74 * qJ(4) + t80 * t75 + pkin(1) + pkin(2);
	t1 = [t72 * t76, t72 * t78, t73, t71 * t79 - t81 * t77 + 0; t73 * t76, t73 * t78, -t72, t71 * t77 + t81 * t79 + 0; -t78, t76, 0, -pkin(4) - qJ(3) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:23:58
	% EndTime: 2020-11-04 21:23:58
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (48->27), mult. (66->31), div. (0->0), fcn. (89->8), ass. (0->19)
	t87 = sin(qJ(6));
	t88 = sin(qJ(5));
	t99 = t87 * t88;
	t90 = cos(qJ(6));
	t98 = t88 * t90;
	t86 = cos(pkin(9));
	t89 = sin(qJ(1));
	t97 = t89 * t86;
	t92 = cos(qJ(1));
	t96 = t92 * t86;
	t85 = sin(pkin(9));
	t93 = pkin(3) + pkin(7);
	t95 = t93 * t85 + qJ(2);
	t91 = cos(qJ(5));
	t94 = pkin(5) * t88 - pkin(8) * t91 + qJ(4);
	t84 = t89 * t85 + t96;
	t83 = -t92 * t85 + t97;
	t82 = t94 * t85 + t93 * t86 + pkin(1) + pkin(2);
	t1 = [-t83 * t98 + t84 * t87, t83 * t99 + t84 * t90, t83 * t91, t82 * t92 + t95 * t89 - t94 * t97 + 0; t83 * t87 + t84 * t98, t83 * t90 - t84 * t99, -t84 * t91, t82 * t89 - t95 * t92 + t94 * t96 + 0; -t91 * t90, t91 * t87, -t88, -t91 * pkin(5) - t88 * pkin(8) - pkin(4) + pkin(6) - qJ(3) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end