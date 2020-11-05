% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRPRP11 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:29
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRRPRP11_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP11_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRP11_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP11_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:29:15
	% EndTime: 2020-11-04 22:29:15
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:29:15
	% EndTime: 2020-11-04 22:29:15
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t90 = cos(qJ(1));
	t89 = sin(qJ(1));
	t1 = [t90, -t89, 0, 0; t89, t90, 0, 0; 0, 0, 1, pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:29:15
	% EndTime: 2020-11-04 22:29:15
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->17), div. (0->0), fcn. (36->6), ass. (0->13)
	t91 = sin(pkin(6));
	t94 = sin(qJ(1));
	t102 = t94 * t91;
	t93 = sin(qJ(2));
	t101 = t94 * t93;
	t95 = cos(qJ(2));
	t100 = t94 * t95;
	t96 = cos(qJ(1));
	t99 = t96 * t91;
	t98 = t96 * t93;
	t97 = t96 * t95;
	t92 = cos(pkin(6));
	t1 = [-t92 * t101 + t97, -t92 * t100 - t98, t102, t96 * pkin(1) + pkin(8) * t102 + 0; t92 * t98 + t100, t92 * t97 - t101, -t99, t94 * pkin(1) - pkin(8) * t99 + 0; t91 * t93, t91 * t95, t92, t92 * pkin(8) + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:29:15
	% EndTime: 2020-11-04 22:29:15
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (29->23), mult. (60->37), div. (0->0), fcn. (81->8), ass. (0->20)
	t106 = sin(pkin(6));
	t111 = cos(qJ(3));
	t121 = t106 * t111;
	t108 = sin(qJ(3));
	t120 = t108 * t106;
	t109 = sin(qJ(2));
	t119 = t109 * t111;
	t110 = sin(qJ(1));
	t118 = t110 * t109;
	t112 = cos(qJ(2));
	t117 = t110 * t112;
	t113 = cos(qJ(1));
	t116 = t113 * t109;
	t115 = t113 * t112;
	t114 = pkin(2) * t109 - pkin(9) * t112;
	t107 = cos(pkin(6));
	t105 = t112 * pkin(2) + t109 * pkin(9) + pkin(1);
	t104 = -t107 * t116 - t117;
	t103 = t106 * pkin(8) - t114 * t107;
	t1 = [(-t107 * t119 + t120) * t110 + t111 * t115, (t107 * t118 - t115) * t108 + t110 * t121, t107 * t117 + t116, t103 * t110 + t105 * t113 + 0; -t104 * t111 - t113 * t120, t104 * t108 - t113 * t121, -t107 * t115 + t118, -t103 * t113 + t105 * t110 + 0; t106 * t119 + t107 * t108, t107 * t111 - t109 * t120, -t106 * t112, t107 * pkin(8) + t114 * t106 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:29:15
	% EndTime: 2020-11-04 22:29:15
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (45->27), mult. (78->41), div. (0->0), fcn. (99->8), ass. (0->22)
	t126 = sin(pkin(6));
	t131 = cos(qJ(3));
	t142 = t126 * t131;
	t128 = sin(qJ(3));
	t141 = t128 * t126;
	t129 = sin(qJ(2));
	t140 = t129 * t131;
	t130 = sin(qJ(1));
	t139 = t130 * t129;
	t132 = cos(qJ(2));
	t138 = t130 * t132;
	t133 = cos(qJ(1));
	t137 = t133 * t129;
	t136 = t133 * t132;
	t125 = t131 * pkin(3) + qJ(4) * t128 + pkin(2);
	t135 = t132 * pkin(9) - t125 * t129;
	t134 = t128 * pkin(3) - qJ(4) * t131 + pkin(8);
	t127 = cos(pkin(6));
	t124 = -t127 * t137 - t138;
	t123 = t129 * pkin(9) + t125 * t132 + pkin(1);
	t122 = t126 * t134 + t135 * t127;
	t1 = [t127 * t138 + t137, (t127 * t140 - t141) * t130 - t131 * t136, (-t127 * t139 + t136) * t128 - t130 * t142, t122 * t130 + t123 * t133 + 0; -t127 * t136 + t139, t124 * t131 + t133 * t141, -t124 * t128 + t133 * t142, -t122 * t133 + t123 * t130 + 0; -t126 * t132, -t126 * t140 - t127 * t128, -t127 * t131 + t129 * t141, -t135 * t126 + t134 * t127 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:29:15
	% EndTime: 2020-11-04 22:29:15
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (70->36), mult. (112->58), div. (0->0), fcn. (146->10), ass. (0->31)
	t150 = cos(pkin(6));
	t156 = cos(qJ(3));
	t172 = t150 * t156;
	t151 = sin(qJ(5));
	t157 = cos(qJ(2));
	t171 = t151 * t157;
	t149 = sin(pkin(6));
	t152 = sin(qJ(3));
	t170 = t152 * t149;
	t153 = sin(qJ(2));
	t169 = t152 * t153;
	t154 = sin(qJ(1));
	t168 = t154 * t157;
	t155 = cos(qJ(5));
	t167 = t155 * t157;
	t166 = t156 * t149;
	t158 = cos(qJ(1));
	t165 = t157 * t158;
	t164 = t158 * t153;
	t160 = pkin(3) + pkin(10);
	t148 = qJ(4) * t152 + t160 * t156 + pkin(2);
	t159 = pkin(4) + pkin(9);
	t163 = t148 * t153 - t159 * t157;
	t162 = qJ(4) * t156 - t160 * t152 - pkin(8);
	t145 = t150 * t169 + t166;
	t161 = -t145 * t151 + t150 * t167;
	t147 = t152 * t171 + t153 * t155;
	t146 = t149 * t169 - t172;
	t144 = t148 * t157 + t159 * t153 + pkin(1);
	t143 = -t149 * t162 - t163 * t150;
	t1 = [t158 * t147 + t161 * t154, (-t145 * t154 + t152 * t165) * t155 - (t150 * t168 + t164) * t151, (-t153 * t172 + t170) * t154 + t156 * t165, t143 * t154 + t144 * t158 + 0; t154 * t147 - t161 * t158, (t145 * t155 + t150 * t171) * t158 + t154 * (-t153 * t151 + t152 * t167), (t150 * t164 + t168) * t156 - t158 * t170, -t143 * t158 + t144 * t154 + 0; t146 * t151 - t149 * t167, t146 * t155 + t149 * t171, t150 * t152 + t153 * t166, t163 * t149 - t162 * t150 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:29:15
	% EndTime: 2020-11-04 22:29:15
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (91->38), mult. (125->60), div. (0->0), fcn. (159->10), ass. (0->32)
	t183 = cos(pkin(6));
	t189 = cos(qJ(3));
	t203 = t183 * t189;
	t184 = sin(qJ(5));
	t190 = cos(qJ(2));
	t202 = t184 * t190;
	t182 = sin(pkin(6));
	t185 = sin(qJ(3));
	t201 = t185 * t182;
	t186 = sin(qJ(2));
	t200 = t185 * t186;
	t187 = sin(qJ(1));
	t199 = t187 * t190;
	t188 = cos(qJ(5));
	t198 = t188 * t190;
	t197 = t189 * t182;
	t191 = cos(qJ(1));
	t196 = t190 * t191;
	t195 = t191 * t186;
	t180 = t184 * pkin(5) + qJ(4);
	t181 = qJ(6) + pkin(3) + pkin(10);
	t175 = t180 * t185 + t181 * t189 + pkin(2);
	t179 = t188 * pkin(5) + pkin(4) + pkin(9);
	t194 = t175 * t186 - t179 * t190;
	t193 = t180 * t189 - t181 * t185 - pkin(8);
	t176 = t183 * t200 + t197;
	t192 = -t176 * t184 + t183 * t198;
	t178 = t185 * t202 + t186 * t188;
	t177 = t182 * t200 - t203;
	t174 = t175 * t190 + t179 * t186 + pkin(1);
	t173 = -t193 * t182 - t194 * t183;
	t1 = [t191 * t178 + t192 * t187, (-t176 * t187 + t185 * t196) * t188 - (t183 * t199 + t195) * t184, (-t186 * t203 + t201) * t187 + t189 * t196, t173 * t187 + t174 * t191 + 0; t187 * t178 - t192 * t191, (t176 * t188 + t183 * t202) * t191 + t187 * (-t186 * t184 + t185 * t198), (t183 * t195 + t199) * t189 - t191 * t201, -t173 * t191 + t174 * t187 + 0; t177 * t184 - t182 * t198, t177 * t188 + t182 * t202, t183 * t185 + t186 * t197, t194 * t182 - t193 * t183 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end