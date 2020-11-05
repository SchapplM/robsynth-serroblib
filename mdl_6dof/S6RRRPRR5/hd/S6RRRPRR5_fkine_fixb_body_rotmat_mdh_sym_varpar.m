% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRPRR5 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:31
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRRPRR5_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR5_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRR5_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR5_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:31:03
	% EndTime: 2020-11-04 22:31:03
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:31:03
	% EndTime: 2020-11-04 22:31:03
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t56 = cos(qJ(1));
	t55 = sin(qJ(1));
	t1 = [t56, -t55, 0, 0; t55, t56, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:31:03
	% EndTime: 2020-11-04 22:31:03
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t60 = cos(qJ(1));
	t59 = cos(qJ(2));
	t58 = sin(qJ(1));
	t57 = sin(qJ(2));
	t1 = [t60 * t59, -t60 * t57, t58, t60 * pkin(1) + t58 * pkin(7) + 0; t58 * t59, -t58 * t57, -t60, t58 * pkin(1) - t60 * pkin(7) + 0; t57, t59, 0, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:31:03
	% EndTime: 2020-11-04 22:31:03
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (19->12), mult. (11->10), div. (0->0), fcn. (19->6), ass. (0->8)
	t67 = pkin(8) + pkin(7);
	t66 = cos(qJ(1));
	t65 = sin(qJ(1));
	t64 = qJ(2) + qJ(3);
	t63 = cos(t64);
	t62 = sin(t64);
	t61 = cos(qJ(2)) * pkin(2) + pkin(1);
	t1 = [t66 * t63, -t66 * t62, t65, t66 * t61 + t67 * t65 + 0; t65 * t63, -t65 * t62, -t66, t65 * t61 - t66 * t67 + 0; t62, t63, 0, sin(qJ(2)) * pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:31:03
	% EndTime: 2020-11-04 22:31:03
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (33->18), mult. (21->14), div. (0->0), fcn. (29->6), ass. (0->8)
	t71 = qJ(2) + qJ(3);
	t69 = sin(t71);
	t70 = cos(t71);
	t75 = pkin(3) * t70 + qJ(4) * t69 + cos(qJ(2)) * pkin(2) + pkin(1);
	t74 = pkin(8) + pkin(7);
	t73 = cos(qJ(1));
	t72 = sin(qJ(1));
	t1 = [t72, -t73 * t70, t73 * t69, t74 * t72 + t75 * t73 + 0; -t73, -t72 * t70, t72 * t69, t75 * t72 - t73 * t74 + 0; 0, -t69, -t70, t69 * pkin(3) - t70 * qJ(4) + sin(qJ(2)) * pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:31:03
	% EndTime: 2020-11-04 22:31:03
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (44->22), mult. (37->25), div. (0->0), fcn. (50->10), ass. (0->18)
	t81 = sin(qJ(5));
	t84 = sin(qJ(1));
	t92 = t84 * t81;
	t85 = cos(qJ(5));
	t91 = t84 * t85;
	t87 = cos(qJ(1));
	t90 = t87 * t81;
	t89 = t87 * t85;
	t88 = pkin(3) + pkin(9);
	t86 = cos(qJ(3));
	t83 = sin(qJ(2));
	t82 = sin(qJ(3));
	t80 = qJ(2) + qJ(3);
	t79 = pkin(4) + pkin(7) + pkin(8);
	t78 = cos(t80);
	t77 = sin(t80);
	t76 = (qJ(4) * t82 + t88 * t86 + pkin(2)) * cos(qJ(2)) + pkin(1) + (qJ(4) * t86 - t82 * t88) * t83;
	t1 = [t77 * t90 + t91, t77 * t89 - t92, t87 * t78, t76 * t87 + t79 * t84 + 0; t77 * t92 - t89, t77 * t91 + t90, t84 * t78, t76 * t84 - t79 * t87 + 0; -t78 * t81, -t78 * t85, t77, t83 * pkin(2) - t78 * qJ(4) + t88 * t77 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:31:03
	% EndTime: 2020-11-04 22:31:04
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (71->30), mult. (51->30), div. (0->0), fcn. (58->14), ass. (0->20)
	t105 = sin(qJ(1));
	t101 = qJ(5) + qJ(6);
	t96 = sin(t101);
	t111 = t105 * t96;
	t98 = cos(t101);
	t110 = t105 * t98;
	t107 = cos(qJ(1));
	t109 = t107 * t96;
	t108 = t107 * t98;
	t102 = qJ(2) + qJ(3);
	t106 = cos(qJ(3));
	t104 = sin(qJ(2));
	t103 = sin(qJ(3));
	t100 = pkin(10) + pkin(3) + pkin(9);
	t99 = cos(t102);
	t97 = sin(t102);
	t95 = sin(qJ(5)) * pkin(5) + qJ(4);
	t94 = cos(qJ(5)) * pkin(5) + pkin(4) + pkin(7) + pkin(8);
	t93 = (t100 * t106 + t95 * t103 + pkin(2)) * cos(qJ(2)) + pkin(1) + (-t103 * t100 + t95 * t106) * t104;
	t1 = [t97 * t109 + t110, t97 * t108 - t111, t107 * t99, t94 * t105 + t93 * t107 + 0; t97 * t111 - t108, t97 * t110 + t109, t105 * t99, t93 * t105 - t94 * t107 + 0; -t99 * t96, -t99 * t98, t97, t100 * t97 - t99 * qJ(4) + t104 * pkin(2) + 0 + pkin(6) + (sin(-qJ(5) + t102) / 0.2e1 - sin(qJ(5) + t102) / 0.2e1) * pkin(5); 0, 0, 0, 1;];
	Tc_mdh = t1;
end