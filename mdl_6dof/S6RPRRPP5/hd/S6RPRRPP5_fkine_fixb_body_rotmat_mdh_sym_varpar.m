% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRRPP5 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:45
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RPRRPP5_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP5_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPP5_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP5_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:45:16
	% EndTime: 2020-11-04 21:45:16
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:45:16
	% EndTime: 2020-11-04 21:45:16
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t63 = cos(qJ(1));
	t62 = sin(qJ(1));
	t1 = [t63, -t62, 0, 0; t62, t63, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:45:16
	% EndTime: 2020-11-04 21:45:16
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t67 = cos(qJ(1));
	t66 = sin(qJ(1));
	t65 = cos(pkin(9));
	t64 = sin(pkin(9));
	t1 = [t67 * t65, -t67 * t64, t66, t67 * pkin(1) + t66 * qJ(2) + 0; t66 * t65, -t66 * t64, -t67, t66 * pkin(1) - t67 * qJ(2) + 0; t64, t65, 0, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:45:16
	% EndTime: 2020-11-04 21:45:16
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (19->12), mult. (11->10), div. (0->0), fcn. (19->6), ass. (0->8)
	t74 = cos(qJ(1));
	t73 = sin(qJ(1));
	t72 = pkin(7) + qJ(2);
	t71 = pkin(9) + qJ(3);
	t70 = cos(t71);
	t69 = sin(t71);
	t68 = cos(pkin(9)) * pkin(2) + pkin(1);
	t1 = [t74 * t70, -t74 * t69, t73, t74 * t68 + t72 * t73 + 0; t73 * t70, -t73 * t69, -t74, t73 * t68 - t74 * t72 + 0; t69, t70, 0, sin(pkin(9)) * pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:45:16
	% EndTime: 2020-11-04 21:45:16
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (37->19), mult. (33->22), div. (0->0), fcn. (46->8), ass. (0->14)
	t80 = sin(qJ(4));
	t81 = sin(qJ(1));
	t88 = t81 * t80;
	t82 = cos(qJ(4));
	t87 = t81 * t82;
	t83 = cos(qJ(1));
	t86 = t83 * t80;
	t85 = t83 * t82;
	t78 = pkin(9) + qJ(3);
	t76 = sin(t78);
	t77 = cos(t78);
	t84 = pkin(3) * t77 + pkin(8) * t76 + cos(pkin(9)) * pkin(2) + pkin(1);
	t79 = pkin(7) + qJ(2);
	t1 = [t77 * t85 + t88, -t77 * t86 + t87, t83 * t76, t79 * t81 + t84 * t83 + 0; t77 * t87 - t86, -t77 * t88 - t85, t81 * t76, -t83 * t79 + t84 * t81 + 0; t76 * t82, -t76 * t80, -t77, t76 * pkin(3) - t77 * pkin(8) + sin(pkin(9)) * pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:45:16
	% EndTime: 2020-11-04 21:45:16
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (52->24), mult. (53->28), div. (0->0), fcn. (70->8), ass. (0->18)
	t98 = sin(qJ(4));
	t99 = sin(qJ(1));
	t106 = t99 * t98;
	t101 = cos(qJ(1));
	t105 = t101 * t98;
	t100 = cos(qJ(4));
	t104 = t99 * t100;
	t103 = t101 * t100;
	t96 = pkin(9) + qJ(3);
	t94 = sin(t96);
	t95 = cos(t96);
	t102 = pkin(3) * t95 + pkin(8) * t94 + cos(pkin(9)) * pkin(2) + pkin(1);
	t97 = pkin(7) + qJ(2);
	t92 = t95 * t103 + t106;
	t91 = t95 * t105 - t104;
	t90 = t95 * t104 - t105;
	t89 = t95 * t106 + t103;
	t1 = [t92, t101 * t94, t91, t92 * pkin(4) + t91 * qJ(5) + t102 * t101 + t97 * t99 + 0; t90, t99 * t94, t89, t90 * pkin(4) + t89 * qJ(5) - t101 * t97 + t102 * t99 + 0; t94 * t100, -t95, t94 * t98, sin(pkin(9)) * pkin(2) - t95 * pkin(8) + pkin(6) + 0 + (pkin(4) * t100 + qJ(5) * t98 + pkin(3)) * t94; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:45:16
	% EndTime: 2020-11-04 21:45:16
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (74->33), mult. (69->36), div. (0->0), fcn. (75->14), ass. (0->21)
	t116 = sin(qJ(4));
	t117 = sin(qJ(1));
	t126 = t117 * t116;
	t118 = cos(qJ(4));
	t125 = t117 * t118;
	t119 = cos(qJ(1));
	t124 = t119 * t116;
	t123 = t119 * t118;
	t112 = pkin(9) + qJ(3);
	t120 = pkin(4) + pkin(5);
	t122 = qJ(5) * t116 + t118 * t120 + pkin(3);
	t121 = qJ(5) * t118 - t120 * t116 - pkin(7) - qJ(2);
	t115 = qJ(6) - pkin(8);
	t114 = cos(pkin(9));
	t113 = sin(pkin(9));
	t111 = -qJ(4) + t112;
	t110 = qJ(4) + t112;
	t109 = cos(t112);
	t108 = sin(t112);
	t107 = (-t113 * t115 + t122 * t114) * cos(qJ(3)) + (-t122 * t113 - t115 * t114) * sin(qJ(3)) + t114 * pkin(2) + pkin(1);
	t1 = [t109 * t123 + t126, t109 * t124 - t125, -t119 * t108, t107 * t119 - t121 * t117 + 0; t109 * t125 - t124, t109 * t126 + t123, -t117 * t108, t107 * t117 + t121 * t119 + 0; t108 * t118, t108 * t116, t109, t115 * t109 + t113 * pkin(2) + t108 * pkin(3) + 0 + pkin(6) + (sin(t111) / 0.2e1 + sin(t110) / 0.2e1) * t120 + (cos(t111) / 0.2e1 - cos(t110) / 0.2e1) * qJ(5); 0, 0, 0, 1;];
	Tc_mdh = t1;
end