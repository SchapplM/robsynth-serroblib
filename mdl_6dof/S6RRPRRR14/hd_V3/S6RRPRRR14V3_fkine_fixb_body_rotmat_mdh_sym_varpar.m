% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPRRR14V3 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:21
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRPRRR14V3_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14V3_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRR14V3_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S6RRPRRR14V3_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [1x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:21:41
	% EndTime: 2020-11-04 22:21:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:21:41
	% EndTime: 2020-11-04 22:21:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (1->1), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t50 = cos(qJ(1));
	t49 = sin(qJ(1));
	t1 = [t50, -t49, 0, 0; t49, t50, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:21:41
	% EndTime: 2020-11-04 22:21:41
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->3), mult. (4->4), div. (0->0), fcn. (12->4), ass. (0->5)
	t54 = cos(qJ(1));
	t53 = cos(qJ(2));
	t52 = sin(qJ(1));
	t51 = sin(qJ(2));
	t1 = [t54 * t53, -t54 * t51, t52, 0; t52 * t53, -t52 * t51, -t54, 0; t51, t53, 0, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:21:41
	% EndTime: 2020-11-04 22:21:41
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (5->5), mult. (9->8), div. (0->0), fcn. (17->4), ass. (0->6)
	t55 = sin(qJ(2));
	t59 = qJ(3) * t55;
	t58 = cos(qJ(1));
	t57 = cos(qJ(2));
	t56 = sin(qJ(1));
	t1 = [t58 * t57, t56, t58 * t55, t58 * t59 + 0; t56 * t57, -t58, t56 * t55, t56 * t59 + 0; t55, 0, -t57, -t57 * qJ(3) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:21:41
	% EndTime: 2020-11-04 22:21:41
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (9->9), mult. (21->17), div. (0->0), fcn. (34->6), ass. (0->11)
	t62 = sin(qJ(1));
	t64 = cos(qJ(2));
	t69 = t62 * t64;
	t60 = sin(qJ(4));
	t65 = cos(qJ(1));
	t68 = t65 * t60;
	t63 = cos(qJ(4));
	t67 = t65 * t63;
	t61 = sin(qJ(2));
	t66 = qJ(3) * t61;
	t1 = [t62 * t60 + t64 * t67, t62 * t63 - t64 * t68, t65 * t61, t65 * t66 + 0; t63 * t69 - t68, -t60 * t69 - t67, t62 * t61, t62 * t66 + 0; t61 * t63, -t61 * t60, -t64, -t64 * qJ(3) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:21:41
	% EndTime: 2020-11-04 22:21:41
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (15->14), mult. (42->31), div. (0->0), fcn. (63->8), ass. (0->17)
	t73 = sin(qJ(2));
	t76 = cos(qJ(4));
	t85 = t73 * t76;
	t78 = cos(qJ(1));
	t84 = t73 * t78;
	t74 = sin(qJ(1));
	t77 = cos(qJ(2));
	t83 = t74 * t77;
	t75 = cos(qJ(5));
	t82 = t77 * t75;
	t72 = sin(qJ(4));
	t81 = t78 * t72;
	t80 = t78 * t76;
	t79 = qJ(3) * t73;
	t71 = sin(qJ(5));
	t70 = t74 * t72 + t77 * t80;
	t1 = [t70 * t75 + t71 * t84, -t70 * t71 + t75 * t84, -t74 * t76 + t77 * t81, t78 * t79 + 0; (t73 * t71 + t76 * t82) * t74 - t75 * t81, (-t76 * t83 + t81) * t71 + t74 * t73 * t75, t72 * t83 + t80, t74 * t79 + 0; -t77 * t71 + t75 * t85, -t71 * t85 - t82, t73 * t72, -t77 * qJ(3) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:21:41
	% EndTime: 2020-11-04 22:21:41
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (28->23), mult. (76->50), div. (0->0), fcn. (110->10), ass. (0->25)
	t91 = sin(qJ(4));
	t95 = cos(qJ(5));
	t109 = t91 * t95;
	t97 = cos(qJ(2));
	t108 = t91 * t97;
	t90 = sin(qJ(5));
	t92 = sin(qJ(2));
	t107 = t92 * t90;
	t106 = t92 * t95;
	t93 = sin(qJ(1));
	t105 = t93 * t97;
	t96 = cos(qJ(4));
	t104 = t95 * t96;
	t103 = t97 * t90;
	t102 = t97 * t95;
	t98 = cos(qJ(1));
	t101 = t98 * t91;
	t100 = t98 * t96;
	t99 = qJ(3) * t92;
	t94 = cos(qJ(6));
	t89 = sin(qJ(6));
	t88 = t96 * t102 + t107;
	t87 = t89 * t109 + t94 * t96;
	t86 = t94 * t108 - t88 * t89;
	t1 = [(t89 * t108 + t88 * t94) * t98 + t93 * (t94 * t109 - t89 * t96), t86 * t98 - t93 * t87, (t97 * t100 + t93 * t91) * t90 - t98 * t106, t98 * t99 + 0; (-t95 * t101 + t88 * t93) * t94 + (t91 * t105 + t100) * t89, t86 * t93 + t98 * t87, (t96 * t105 - t101) * t90 - t93 * t106, t93 * t99 + 0; (t94 * t104 + t89 * t91) * t92 - t94 * t103, (-t92 * t104 + t103) * t89 + t92 * t91 * t94, t96 * t107 + t102, -t97 * qJ(3) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end