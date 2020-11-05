% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPRRR12 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:20
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRPRRR12_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR12_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRR12_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR12_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:20:26
	% EndTime: 2020-11-04 22:20:26
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:20:26
	% EndTime: 2020-11-04 22:20:26
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t77 = cos(qJ(1));
	t76 = sin(qJ(1));
	t1 = [t77, -t76, 0, 0; t76, t77, 0, 0; 0, 0, 1, pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:20:26
	% EndTime: 2020-11-04 22:20:26
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->17), div. (0->0), fcn. (36->6), ass. (0->13)
	t78 = sin(pkin(6));
	t81 = sin(qJ(1));
	t89 = t81 * t78;
	t80 = sin(qJ(2));
	t88 = t81 * t80;
	t82 = cos(qJ(2));
	t87 = t81 * t82;
	t83 = cos(qJ(1));
	t86 = t83 * t78;
	t85 = t83 * t80;
	t84 = t83 * t82;
	t79 = cos(pkin(6));
	t1 = [-t79 * t88 + t84, -t79 * t87 - t85, t89, t83 * pkin(1) + pkin(8) * t89 + 0; t79 * t85 + t87, t79 * t84 - t88, -t86, t81 * pkin(1) - pkin(8) * t86 + 0; t78 * t80, t78 * t82, t79, t79 * pkin(8) + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:20:26
	% EndTime: 2020-11-04 22:20:26
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (23->18), mult. (39->24), div. (0->0), fcn. (52->6), ass. (0->14)
	t94 = sin(qJ(2));
	t95 = sin(qJ(1));
	t102 = t95 * t94;
	t96 = cos(qJ(2));
	t101 = t95 * t96;
	t97 = cos(qJ(1));
	t100 = t97 * t94;
	t99 = t97 * t96;
	t98 = pkin(2) * t94 - qJ(3) * t96;
	t93 = cos(pkin(6));
	t92 = sin(pkin(6));
	t91 = t96 * pkin(2) + t94 * qJ(3) + pkin(1);
	t90 = t92 * pkin(8) - t98 * t93;
	t1 = [t95 * t92, t93 * t102 - t99, t93 * t101 + t100, t90 * t95 + t91 * t97 + 0; -t97 * t92, -t93 * t100 - t101, -t93 * t99 + t102, -t90 * t97 + t91 * t95 + 0; t93, -t92 * t94, -t92 * t96, t93 * pkin(8) + t98 * t92 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:20:26
	% EndTime: 2020-11-04 22:20:26
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (36->24), mult. (60->37), div. (0->0), fcn. (81->8), ass. (0->22)
	t106 = sin(pkin(6));
	t111 = cos(qJ(4));
	t123 = t106 * t111;
	t108 = sin(qJ(4));
	t122 = t108 * t106;
	t109 = sin(qJ(2));
	t110 = sin(qJ(1));
	t121 = t110 * t109;
	t112 = cos(qJ(2));
	t120 = t110 * t112;
	t119 = t111 * t112;
	t113 = cos(qJ(1));
	t118 = t113 * t109;
	t117 = t113 * t112;
	t115 = pkin(2) + pkin(9);
	t116 = qJ(3) * t112 - t109 * t115;
	t114 = pkin(3) + pkin(8);
	t107 = cos(pkin(6));
	t105 = t109 * qJ(3) + t115 * t112 + pkin(1);
	t104 = -t107 * t117 + t121;
	t103 = t106 * t114 + t116 * t107;
	t1 = [t110 * t123 + (t107 * t120 + t118) * t108, (t107 * t119 - t122) * t110 + t111 * t118, -t107 * t121 + t117, t103 * t110 + t105 * t113 + 0; t104 * t108 - t113 * t123, t104 * t111 + t113 * t122, t107 * t118 + t120, -t103 * t113 + t105 * t110 + 0; t107 * t111 - t112 * t122, -t106 * t119 - t107 * t108, t106 * t109, -t116 * t106 + t114 * t107 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:20:26
	% EndTime: 2020-11-04 22:20:26
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (61->26), mult. (69->38), div. (0->0), fcn. (90->10), ass. (0->25)
	t133 = sin(pkin(6));
	t136 = sin(qJ(1));
	t148 = t133 * t136;
	t138 = cos(qJ(2));
	t147 = t133 * t138;
	t139 = cos(qJ(1));
	t146 = t133 * t139;
	t135 = sin(qJ(2));
	t145 = t136 * t135;
	t144 = t136 * t138;
	t143 = t139 * t135;
	t142 = t139 * t138;
	t141 = cos(qJ(4)) * pkin(4) + pkin(3) + pkin(8);
	t128 = sin(qJ(4)) * pkin(4) + qJ(3);
	t131 = pkin(2) + pkin(9) + pkin(10);
	t140 = t128 * t138 - t131 * t135;
	t134 = cos(pkin(6));
	t132 = qJ(4) + qJ(5);
	t130 = cos(t132);
	t129 = sin(t132);
	t127 = t134 * t144 + t143;
	t126 = -t134 * t142 + t145;
	t125 = t128 * t135 + t131 * t138 + pkin(1);
	t124 = t133 * t141 + t140 * t134;
	t1 = [t127 * t129 + t130 * t148, t127 * t130 - t129 * t148, -t134 * t145 + t142, t124 * t136 + t125 * t139 + 0; t126 * t129 - t130 * t146, t126 * t130 + t129 * t146, t134 * t143 + t144, -t124 * t139 + t125 * t136 + 0; -t129 * t147 + t134 * t130, -t134 * t129 - t130 * t147, t133 * t135, -t140 * t133 + t141 * t134 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:20:26
	% EndTime: 2020-11-04 22:20:26
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (108->39), mult. (139->56), div. (0->0), fcn. (183->12), ass. (0->36)
	t166 = sin(pkin(6));
	t169 = sin(qJ(2));
	t184 = t166 * t169;
	t170 = sin(qJ(1));
	t183 = t166 * t170;
	t173 = cos(qJ(2));
	t182 = t166 * t173;
	t174 = cos(qJ(1));
	t181 = t166 * t174;
	t180 = t170 * t169;
	t179 = t170 * t173;
	t178 = t174 * t169;
	t177 = t174 * t173;
	t176 = cos(qJ(4)) * pkin(4) + pkin(3) + pkin(8);
	t161 = sin(qJ(4)) * pkin(4) + qJ(3);
	t164 = pkin(2) + pkin(9) + pkin(10);
	t175 = t161 * t173 - t164 * t169;
	t171 = cos(qJ(6));
	t168 = sin(qJ(6));
	t167 = cos(pkin(6));
	t165 = qJ(4) + qJ(5);
	t163 = cos(t165);
	t162 = sin(t165);
	t160 = -t167 * t180 + t177;
	t159 = t167 * t179 + t178;
	t158 = t167 * t178 + t179;
	t157 = -t167 * t177 + t180;
	t156 = t161 * t169 + t164 * t173 + pkin(1);
	t155 = -t162 * t182 + t163 * t167;
	t154 = t162 * t167 + t163 * t182;
	t153 = t157 * t162 - t163 * t181;
	t152 = t157 * t163 + t162 * t181;
	t151 = t159 * t162 + t163 * t183;
	t150 = -t159 * t163 + t162 * t183;
	t149 = t166 * t176 + t167 * t175;
	t1 = [t151 * t171 + t160 * t168, -t151 * t168 + t160 * t171, t150, pkin(5) * t151 + pkin(11) * t150 + t149 * t170 + t156 * t174 + 0; t153 * t171 + t158 * t168, -t153 * t168 + t158 * t171, -t152, pkin(5) * t153 - pkin(11) * t152 - t149 * t174 + t156 * t170 + 0; t155 * t171 + t168 * t184, -t155 * t168 + t171 * t184, t154, t155 * pkin(5) + t154 * pkin(11) - t166 * t175 + t167 * t176 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end