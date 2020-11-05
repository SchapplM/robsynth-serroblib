% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRPRP4 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:27
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRRPRP4_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP4_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRP4_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP4_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:27:05
	% EndTime: 2020-11-04 22:27:05
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:27:05
	% EndTime: 2020-11-04 22:27:05
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t59 = cos(qJ(1));
	t58 = sin(qJ(1));
	t1 = [t59, -t58, 0, 0; t58, t59, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:27:05
	% EndTime: 2020-11-04 22:27:05
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t63 = cos(qJ(1));
	t62 = cos(qJ(2));
	t61 = sin(qJ(1));
	t60 = sin(qJ(2));
	t1 = [t63 * t62, -t63 * t60, t61, t63 * pkin(1) + t61 * pkin(7) + 0; t61 * t62, -t61 * t60, -t63, t61 * pkin(1) - t63 * pkin(7) + 0; t60, t62, 0, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:27:05
	% EndTime: 2020-11-04 22:27:05
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (19->12), mult. (11->10), div. (0->0), fcn. (19->6), ass. (0->8)
	t70 = pkin(8) + pkin(7);
	t69 = cos(qJ(1));
	t68 = sin(qJ(1));
	t67 = qJ(2) + qJ(3);
	t66 = cos(t67);
	t65 = sin(t67);
	t64 = cos(qJ(2)) * pkin(2) + pkin(1);
	t1 = [t69 * t66, -t69 * t65, t68, t69 * t64 + t70 * t68 + 0; t68 * t66, -t68 * t65, -t69, t68 * t64 - t69 * t70 + 0; t65, t66, 0, sin(qJ(2)) * pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:27:05
	% EndTime: 2020-11-04 22:27:05
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (33->18), mult. (21->14), div. (0->0), fcn. (29->6), ass. (0->8)
	t74 = qJ(2) + qJ(3);
	t72 = sin(t74);
	t73 = cos(t74);
	t78 = pkin(3) * t73 + qJ(4) * t72 + cos(qJ(2)) * pkin(2) + pkin(1);
	t77 = pkin(8) + pkin(7);
	t76 = cos(qJ(1));
	t75 = sin(qJ(1));
	t1 = [t75, -t76 * t73, t76 * t72, t77 * t75 + t76 * t78 + 0; -t76, -t75 * t73, t75 * t72, t75 * t78 - t76 * t77 + 0; 0, -t72, -t73, t72 * pkin(3) - t73 * qJ(4) + sin(qJ(2)) * pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:27:05
	% EndTime: 2020-11-04 22:27:05
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (44->22), mult. (37->25), div. (0->0), fcn. (50->10), ass. (0->18)
	t84 = sin(qJ(5));
	t87 = sin(qJ(1));
	t95 = t87 * t84;
	t88 = cos(qJ(5));
	t94 = t87 * t88;
	t90 = cos(qJ(1));
	t93 = t90 * t84;
	t92 = t90 * t88;
	t91 = pkin(3) + pkin(9);
	t89 = cos(qJ(3));
	t86 = sin(qJ(2));
	t85 = sin(qJ(3));
	t83 = qJ(2) + qJ(3);
	t82 = pkin(4) + pkin(7) + pkin(8);
	t81 = cos(t83);
	t80 = sin(t83);
	t79 = (qJ(4) * t85 + t91 * t89 + pkin(2)) * cos(qJ(2)) + pkin(1) + (qJ(4) * t89 - t85 * t91) * t86;
	t1 = [t80 * t93 + t94, t80 * t92 - t95, t90 * t81, t79 * t90 + t82 * t87 + 0; t80 * t95 - t92, t80 * t94 + t93, t87 * t81, t79 * t87 - t82 * t90 + 0; -t81 * t84, -t81 * t88, t80, t86 * pkin(2) - t81 * qJ(4) + t91 * t80 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:27:05
	% EndTime: 2020-11-04 22:27:05
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (59->28), mult. (57->31), div. (0->0), fcn. (74->10), ass. (0->22)
	t105 = sin(qJ(5));
	t108 = sin(qJ(1));
	t116 = t108 * t105;
	t109 = cos(qJ(5));
	t115 = t108 * t109;
	t111 = cos(qJ(1));
	t114 = t111 * t105;
	t113 = t111 * t109;
	t112 = pkin(3) + pkin(9);
	t110 = cos(qJ(3));
	t107 = sin(qJ(2));
	t106 = sin(qJ(3));
	t104 = qJ(2) + qJ(3);
	t103 = pkin(4) + pkin(7) + pkin(8);
	t102 = cos(t104);
	t101 = sin(t104);
	t100 = t101 * t116 - t113;
	t99 = t101 * t115 + t114;
	t98 = t101 * t114 + t115;
	t97 = -t101 * t113 + t116;
	t96 = (qJ(4) * t106 + t112 * t110 + pkin(2)) * cos(qJ(2)) + pkin(1) + (qJ(4) * t110 - t106 * t112) * t107;
	t1 = [t98, t111 * t102, t97, t98 * pkin(5) + t97 * qJ(6) + t103 * t108 + t96 * t111 + 0; t100, t108 * t102, -t99, t100 * pkin(5) - t99 * qJ(6) - t103 * t111 + t96 * t108 + 0; -t102 * t105, t101, t102 * t109, t107 * pkin(2) + t112 * t101 + pkin(6) + 0 + (-pkin(5) * t105 + qJ(6) * t109 - qJ(4)) * t102; 0, 0, 0, 1;];
	Tc_mdh = t1;
end