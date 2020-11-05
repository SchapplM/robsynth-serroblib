% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRRRRR3 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:21
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6PRRRRR3_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR3_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRRR3_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR3_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [12x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:21:09
	% EndTime: 2020-11-04 21:21:09
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:21:09
	% EndTime: 2020-11-04 21:21:09
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t98 = cos(pkin(12));
	t97 = sin(pkin(12));
	t1 = [t98, -t97, 0, 0; t97, t98, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:21:09
	% EndTime: 2020-11-04 21:21:09
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->19), div. (0->0), fcn. (36->6), ass. (0->11)
	t100 = sin(pkin(6));
	t99 = sin(pkin(12));
	t108 = t99 * t100;
	t101 = cos(pkin(12));
	t107 = t101 * t100;
	t102 = cos(pkin(6));
	t103 = sin(qJ(2));
	t106 = t102 * t103;
	t104 = cos(qJ(2));
	t105 = t102 * t104;
	t1 = [t101 * t104 - t99 * t106, -t101 * t103 - t99 * t105, t108, t101 * pkin(1) + pkin(7) * t108 + 0; t101 * t106 + t99 * t104, t101 * t105 - t99 * t103, -t107, t99 * pkin(1) - pkin(7) * t107 + 0; t100 * t103, t100 * t104, t102, t102 * pkin(7) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:21:09
	% EndTime: 2020-11-04 21:21:09
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (29->27), mult. (64->48), div. (0->0), fcn. (85->8), ass. (0->18)
	t112 = sin(pkin(6));
	t125 = t112 * pkin(7);
	t111 = sin(pkin(12));
	t114 = cos(pkin(6));
	t124 = t111 * t114;
	t115 = sin(qJ(3));
	t123 = t112 * t115;
	t117 = cos(qJ(3));
	t122 = t112 * t117;
	t113 = cos(pkin(12));
	t121 = t113 * t114;
	t116 = sin(qJ(2));
	t120 = t114 * t116;
	t118 = cos(qJ(2));
	t119 = t114 * t118;
	t110 = -t111 * t120 + t113 * t118;
	t109 = t111 * t118 + t113 * t120;
	t1 = [t110 * t117 + t111 * t123, -t110 * t115 + t111 * t122, t111 * t119 + t113 * t116, (t113 * pkin(2) + pkin(8) * t124) * t118 + (-pkin(2) * t124 + t113 * pkin(8)) * t116 + t111 * t125 + t113 * pkin(1) + 0; t109 * t117 - t113 * t123, -t109 * t115 - t113 * t122, t111 * t116 - t113 * t119, (t111 * pkin(2) - pkin(8) * t121) * t118 + (pkin(2) * t121 + t111 * pkin(8)) * t116 - t113 * t125 + t111 * pkin(1) + 0; t114 * t115 + t116 * t122, t114 * t117 - t116 * t123, -t112 * t118, t114 * pkin(7) + qJ(1) + 0 + (pkin(2) * t116 - pkin(8) * t118) * t112; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:21:09
	% EndTime: 2020-11-04 21:21:09
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (57->32), mult. (136->56), div. (0->0), fcn. (170->10), ass. (0->29)
	t133 = sin(pkin(6));
	t135 = cos(pkin(6));
	t138 = sin(qJ(2));
	t141 = cos(qJ(2));
	t137 = sin(qJ(3));
	t140 = cos(qJ(3));
	t145 = t140 * pkin(3) + t137 * pkin(9) + pkin(2);
	t143 = -pkin(8) * t141 + t145 * t138;
	t144 = t137 * pkin(3) - t140 * pkin(9) + pkin(7);
	t156 = t144 * t133 - t143 * t135;
	t152 = t133 * t140;
	t151 = t133 * t141;
	t150 = t135 * t138;
	t149 = t135 * t141;
	t148 = t137 * t133;
	t147 = t138 * t140;
	t146 = t140 * t141;
	t142 = pkin(8) * t138 + t145 * t141 + pkin(1);
	t139 = cos(qJ(4));
	t136 = sin(qJ(4));
	t134 = cos(pkin(12));
	t132 = sin(pkin(12));
	t131 = t135 * t147 - t148;
	t130 = t133 * t147 + t135 * t137;
	t129 = t132 * t149 + t134 * t138;
	t128 = t132 * t138 - t134 * t149;
	t127 = -t132 * t131 + t134 * t146;
	t126 = t134 * t131 + t132 * t146;
	t1 = [t127 * t139 + t129 * t136, -t127 * t136 + t129 * t139, -(t132 * t150 - t134 * t141) * t137 - t132 * t152, t156 * t132 + t142 * t134 + 0; t126 * t139 + t128 * t136, -t126 * t136 + t128 * t139, (t132 * t141 + t134 * t150) * t137 + t134 * t152, t142 * t132 - t156 * t134 + 0; t130 * t139 - t136 * t151, -t130 * t136 - t139 * t151, -t135 * t140 + t138 * t148, t143 * t133 + t144 * t135 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:21:09
	% EndTime: 2020-11-04 21:21:09
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (80->45), mult. (150->70), div. (0->0), fcn. (196->12), ass. (0->33)
	t189 = pkin(4) * sin(qJ(4));
	t172 = sin(pkin(6));
	t188 = t172 * pkin(7);
	t171 = sin(pkin(12));
	t174 = cos(pkin(6));
	t187 = t171 * t174;
	t176 = sin(qJ(3));
	t186 = t172 * t176;
	t178 = cos(qJ(3));
	t185 = t172 * t178;
	t179 = cos(qJ(2));
	t184 = t172 * t179;
	t173 = cos(pkin(12));
	t183 = t173 * t174;
	t177 = sin(qJ(2));
	t182 = t174 * t177;
	t181 = t174 * t179;
	t180 = -pkin(10) - pkin(9);
	t170 = qJ(4) + qJ(5);
	t169 = cos(t170);
	t168 = sin(t170);
	t167 = cos(qJ(4)) * pkin(4) + pkin(3);
	t166 = t174 * t176 + t177 * t185;
	t165 = -t174 * t178 + t177 * t186;
	t164 = -t171 * t182 + t173 * t179;
	t163 = t171 * t181 + t173 * t177;
	t162 = t171 * t179 + t173 * t182;
	t161 = t171 * t177 - t173 * t181;
	t160 = -t164 * t176 + t171 * t185;
	t159 = t164 * t178 + t171 * t186;
	t158 = t162 * t178 - t173 * t186;
	t157 = t162 * t176 + t173 * t185;
	t1 = [t159 * t169 + t163 * t168, -t159 * t168 + t163 * t169, -t160, t159 * t167 + t160 * t180 + t163 * t189 + (t173 * pkin(2) + pkin(8) * t187) * t179 + (-pkin(2) * t187 + t173 * pkin(8)) * t177 + t171 * t188 + t173 * pkin(1) + 0; t158 * t169 + t161 * t168, -t158 * t168 + t161 * t169, t157, t158 * t167 - t157 * t180 + t161 * t189 + (t171 * pkin(2) - pkin(8) * t183) * t179 + (pkin(2) * t183 + t171 * pkin(8)) * t177 - t173 * t188 + t171 * pkin(1) + 0; t166 * t169 - t168 * t184, -t166 * t168 - t169 * t184, t165, t174 * pkin(7) - t165 * t180 + t166 * t167 + qJ(1) + 0 + (pkin(2) * t177 + (-pkin(8) - t189) * t179) * t172; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:21:09
	% EndTime: 2020-11-04 21:21:09
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (107->48), mult. (156->72), div. (0->0), fcn. (202->14), ass. (0->34)
	t208 = sin(pkin(6));
	t222 = t208 * pkin(7);
	t207 = sin(pkin(12));
	t210 = cos(pkin(6));
	t221 = t207 * t210;
	t211 = sin(qJ(3));
	t220 = t208 * t211;
	t213 = cos(qJ(3));
	t219 = t208 * t213;
	t214 = cos(qJ(2));
	t218 = t208 * t214;
	t209 = cos(pkin(12));
	t217 = t209 * t210;
	t212 = sin(qJ(2));
	t216 = t210 * t212;
	t215 = t210 * t214;
	t206 = qJ(4) + qJ(5);
	t205 = -pkin(11) - pkin(10) - pkin(9);
	t204 = qJ(6) + t206;
	t203 = cos(t204);
	t202 = sin(t204);
	t201 = pkin(5) * sin(t206) + pkin(4) * sin(qJ(4));
	t200 = pkin(5) * cos(t206) + cos(qJ(4)) * pkin(4) + pkin(3);
	t199 = t210 * t211 + t212 * t219;
	t198 = -t210 * t213 + t212 * t220;
	t197 = -t207 * t216 + t209 * t214;
	t196 = t207 * t215 + t209 * t212;
	t195 = t207 * t214 + t209 * t216;
	t194 = t207 * t212 - t209 * t215;
	t193 = -t197 * t211 + t207 * t219;
	t192 = t197 * t213 + t207 * t220;
	t191 = t195 * t213 - t209 * t220;
	t190 = t195 * t211 + t209 * t219;
	t1 = [t192 * t203 + t196 * t202, -t192 * t202 + t196 * t203, -t193, t192 * t200 + t193 * t205 + t196 * t201 + (t209 * pkin(2) + pkin(8) * t221) * t214 + (-pkin(2) * t221 + t209 * pkin(8)) * t212 + t207 * t222 + t209 * pkin(1) + 0; t191 * t203 + t194 * t202, -t191 * t202 + t194 * t203, t190, t191 * t200 - t190 * t205 + t194 * t201 + (t207 * pkin(2) - pkin(8) * t217) * t214 + (pkin(2) * t217 + t207 * pkin(8)) * t212 - t209 * t222 + t207 * pkin(1) + 0; t199 * t203 - t202 * t218, -t199 * t202 - t203 * t218, t198, t210 * pkin(7) - t198 * t205 + t199 * t200 + qJ(1) + 0 + (pkin(2) * t212 + (-pkin(8) - t201) * t214) * t208; 0, 0, 0, 1;];
	Tc_mdh = t1;
end