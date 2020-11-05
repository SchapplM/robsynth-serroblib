% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRRPR6 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:48
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RPRRPR6_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR6_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPR6_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR6_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:48:19
	% EndTime: 2020-11-04 21:48:19
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:48:19
	% EndTime: 2020-11-04 21:48:19
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t59 = cos(qJ(1));
	t58 = sin(qJ(1));
	t1 = [t59, -t58, 0, 0; t58, t59, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:48:19
	% EndTime: 2020-11-04 21:48:19
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t63 = cos(qJ(1));
	t62 = sin(qJ(1));
	t61 = cos(pkin(10));
	t60 = sin(pkin(10));
	t1 = [t63 * t61, -t63 * t60, t62, t63 * pkin(1) + t62 * qJ(2) + 0; t62 * t61, -t62 * t60, -t63, t62 * pkin(1) - t63 * qJ(2) + 0; t60, t61, 0, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:48:19
	% EndTime: 2020-11-04 21:48:20
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (19->12), mult. (11->10), div. (0->0), fcn. (19->6), ass. (0->8)
	t70 = cos(qJ(1));
	t69 = sin(qJ(1));
	t68 = pkin(7) + qJ(2);
	t67 = pkin(10) + qJ(3);
	t66 = cos(t67);
	t65 = sin(t67);
	t64 = cos(pkin(10)) * pkin(2) + pkin(1);
	t1 = [t70 * t66, -t70 * t65, t69, t70 * t64 + t68 * t69 + 0; t69 * t66, -t69 * t65, -t70, t69 * t64 - t70 * t68 + 0; t65, t66, 0, sin(pkin(10)) * pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:48:20
	% EndTime: 2020-11-04 21:48:20
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (37->19), mult. (33->22), div. (0->0), fcn. (46->8), ass. (0->14)
	t76 = sin(qJ(4));
	t77 = sin(qJ(1));
	t84 = t77 * t76;
	t78 = cos(qJ(4));
	t83 = t77 * t78;
	t79 = cos(qJ(1));
	t82 = t79 * t76;
	t81 = t79 * t78;
	t74 = pkin(10) + qJ(3);
	t72 = sin(t74);
	t73 = cos(t74);
	t80 = pkin(3) * t73 + pkin(8) * t72 + cos(pkin(10)) * pkin(2) + pkin(1);
	t75 = pkin(7) + qJ(2);
	t1 = [t73 * t81 + t84, -t73 * t82 + t83, t79 * t72, t75 * t77 + t80 * t79 + 0; t73 * t83 - t82, -t73 * t84 - t81, t77 * t72, -t79 * t75 + t80 * t77 + 0; t72 * t78, -t72 * t76, -t73, t72 * pkin(3) - t73 * pkin(8) + sin(pkin(10)) * pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:48:20
	% EndTime: 2020-11-04 21:48:20
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (55->23), mult. (40->24), div. (0->0), fcn. (53->10), ass. (0->17)
	t92 = qJ(4) + pkin(11);
	t88 = sin(t92);
	t96 = sin(qJ(1));
	t103 = t96 * t88;
	t90 = cos(t92);
	t102 = t96 * t90;
	t97 = cos(qJ(1));
	t101 = t97 * t88;
	t100 = t97 * t90;
	t99 = pkin(4) * sin(qJ(4)) + pkin(7) + qJ(2);
	t86 = cos(qJ(4)) * pkin(4) + pkin(3);
	t91 = pkin(10) + qJ(3);
	t87 = sin(t91);
	t89 = cos(t91);
	t93 = -qJ(5) - pkin(8);
	t98 = t86 * t89 - t87 * t93 + cos(pkin(10)) * pkin(2) + pkin(1);
	t1 = [t89 * t100 + t103, -t89 * t101 + t102, t97 * t87, t99 * t96 + t98 * t97 + 0; t89 * t102 - t101, -t89 * t103 - t100, t96 * t87, t98 * t96 - t99 * t97 + 0; t87 * t90, -t87 * t88, -t89, t87 * t86 + t89 * t93 + sin(pkin(10)) * pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:48:20
	% EndTime: 2020-11-04 21:48:20
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (78->26), mult. (45->26), div. (0->0), fcn. (58->12), ass. (0->18)
	t114 = qJ(4) + pkin(11);
	t111 = qJ(6) + t114;
	t106 = sin(t111);
	t116 = sin(qJ(1));
	t123 = t116 * t106;
	t107 = cos(t111);
	t122 = t116 * t107;
	t117 = cos(qJ(1));
	t121 = t117 * t106;
	t120 = t117 * t107;
	t119 = pkin(5) * sin(t114) + pkin(4) * sin(qJ(4)) + pkin(7) + qJ(2);
	t104 = pkin(5) * cos(t114) + cos(qJ(4)) * pkin(4) + pkin(3);
	t113 = pkin(10) + qJ(3);
	t109 = sin(t113);
	t110 = cos(t113);
	t112 = -pkin(9) - qJ(5) - pkin(8);
	t118 = t104 * t110 - t109 * t112 + cos(pkin(10)) * pkin(2) + pkin(1);
	t1 = [t110 * t120 + t123, -t110 * t121 + t122, t117 * t109, t119 * t116 + t118 * t117 + 0; t110 * t122 - t121, -t110 * t123 - t120, t116 * t109, t118 * t116 - t119 * t117 + 0; t109 * t107, -t109 * t106, -t110, t109 * t104 + t110 * t112 + sin(pkin(10)) * pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end