% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRPRR14 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5,theta3]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:38
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RRPRR14_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR14_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPRR14_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRPRR14_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:38:47
	% EndTime: 2020-11-04 20:38:47
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:38:47
	% EndTime: 2020-11-04 20:38:47
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t68 = cos(qJ(1));
	t67 = sin(qJ(1));
	t1 = [t68, -t67, 0, 0; t67, t68, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:38:47
	% EndTime: 2020-11-04 20:38:47
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (11->11), mult. (23->17), div. (0->0), fcn. (36->6), ass. (0->13)
	t69 = sin(pkin(5));
	t72 = sin(qJ(1));
	t80 = t72 * t69;
	t71 = sin(qJ(2));
	t79 = t72 * t71;
	t73 = cos(qJ(2));
	t78 = t72 * t73;
	t74 = cos(qJ(1));
	t77 = t74 * t69;
	t76 = t74 * t71;
	t75 = t74 * t73;
	t70 = cos(pkin(5));
	t1 = [-t70 * t79 + t75, -t70 * t78 - t76, t80, t74 * pkin(1) + pkin(7) * t80 + 0; t70 * t76 + t78, t70 * t75 - t79, -t77, t72 * pkin(1) - pkin(7) * t77 + 0; t69 * t71, t69 * t73, t70, t70 * pkin(7) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:38:47
	% EndTime: 2020-11-04 20:38:47
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (29->22), mult. (60->36), div. (0->0), fcn. (81->8), ass. (0->21)
	t86 = sin(pkin(5));
	t89 = sin(qJ(2));
	t100 = t86 * t89;
	t90 = sin(qJ(1));
	t99 = t86 * t90;
	t92 = cos(qJ(1));
	t98 = t86 * t92;
	t97 = t90 * t89;
	t91 = cos(qJ(2));
	t96 = t90 * t91;
	t95 = t92 * t89;
	t94 = t92 * t91;
	t93 = pkin(2) * t89 - qJ(3) * t91;
	t88 = cos(pkin(5));
	t87 = cos(pkin(10));
	t85 = sin(pkin(10));
	t84 = t91 * pkin(2) + t89 * qJ(3) + pkin(1);
	t83 = t88 * t95 + t96;
	t82 = t88 * t97 - t94;
	t81 = t86 * pkin(7) - t93 * t88;
	t1 = [-t82 * t87 + t85 * t99, t82 * t85 + t87 * t99, t88 * t96 + t95, t81 * t90 + t84 * t92 + 0; t83 * t87 - t85 * t98, -t83 * t85 - t87 * t98, -t88 * t94 + t97, -t81 * t92 + t84 * t90 + 0; t87 * t100 + t88 * t85, -t85 * t100 + t88 * t87, -t86 * t91, t88 * pkin(7) + t93 * t86 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:38:47
	% EndTime: 2020-11-04 20:38:48
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (54->26), mult. (68->38), div. (0->0), fcn. (89->10), ass. (0->25)
	t110 = sin(pkin(5));
	t113 = sin(qJ(2));
	t124 = t110 * t113;
	t114 = sin(qJ(1));
	t123 = t110 * t114;
	t116 = cos(qJ(1));
	t122 = t110 * t116;
	t121 = t114 * t113;
	t115 = cos(qJ(2));
	t120 = t114 * t115;
	t119 = t116 * t113;
	t118 = t116 * t115;
	t106 = cos(pkin(10)) * pkin(3) + pkin(2);
	t112 = qJ(3) + pkin(8);
	t117 = t106 * t113 - t112 * t115;
	t111 = cos(pkin(5));
	t109 = pkin(10) + qJ(4);
	t108 = cos(t109);
	t107 = sin(t109);
	t105 = sin(pkin(10)) * pkin(3) + pkin(7);
	t104 = t111 * t119 + t120;
	t103 = t111 * t121 - t118;
	t102 = t106 * t115 + t112 * t113 + pkin(1);
	t101 = t110 * t105 - t117 * t111;
	t1 = [-t103 * t108 + t107 * t123, t103 * t107 + t108 * t123, t111 * t120 + t119, t101 * t114 + t102 * t116 + 0; t104 * t108 - t107 * t122, -t104 * t107 - t108 * t122, -t111 * t118 + t121, -t101 * t116 + t102 * t114 + 0; t111 * t107 + t108 * t124, -t107 * t124 + t111 * t108, -t110 * t115, t105 * t111 + t117 * t110 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:38:48
	% EndTime: 2020-11-04 20:38:48
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (100->38), mult. (138->56), div. (0->0), fcn. (182->12), ass. (0->36)
	t142 = sin(pkin(5));
	t146 = sin(qJ(2));
	t159 = t142 * t146;
	t147 = sin(qJ(1));
	t158 = t142 * t147;
	t149 = cos(qJ(2));
	t157 = t142 * t149;
	t150 = cos(qJ(1));
	t156 = t142 * t150;
	t155 = t147 * t146;
	t154 = t147 * t149;
	t153 = t150 * t146;
	t152 = t150 * t149;
	t138 = cos(pkin(10)) * pkin(3) + pkin(2);
	t144 = qJ(3) + pkin(8);
	t151 = t138 * t146 - t144 * t149;
	t148 = cos(qJ(5));
	t145 = sin(qJ(5));
	t143 = cos(pkin(5));
	t141 = pkin(10) + qJ(4);
	t140 = cos(t141);
	t139 = sin(t141);
	t137 = sin(pkin(10)) * pkin(3) + pkin(7);
	t136 = t143 * t154 + t153;
	t135 = t143 * t153 + t154;
	t134 = -t143 * t152 + t155;
	t133 = t143 * t155 - t152;
	t132 = t138 * t149 + t144 * t146 + pkin(1);
	t131 = t143 * t139 + t140 * t159;
	t130 = t139 * t159 - t143 * t140;
	t129 = t142 * t137 - t151 * t143;
	t128 = -t133 * t140 + t139 * t158;
	t127 = t135 * t140 - t139 * t156;
	t126 = t135 * t139 + t140 * t156;
	t125 = t133 * t139 + t140 * t158;
	t1 = [t128 * t148 + t136 * t145, -t128 * t145 + t136 * t148, -t125, t128 * pkin(4) - t125 * pkin(9) + t129 * t147 + t132 * t150 + 0; t127 * t148 + t134 * t145, -t127 * t145 + t134 * t148, t126, t127 * pkin(4) + t126 * pkin(9) - t129 * t150 + t132 * t147 + 0; t131 * t148 - t145 * t157, -t131 * t145 - t148 * t157, t130, t131 * pkin(4) + t130 * pkin(9) + t137 * t143 + t151 * t142 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end