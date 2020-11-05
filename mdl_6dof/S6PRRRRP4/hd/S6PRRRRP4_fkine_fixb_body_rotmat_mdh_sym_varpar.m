% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRRRRP4 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:19
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6PRRRRP4_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP4_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRRP4_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP4_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:19:02
	% EndTime: 2020-11-04 21:19:02
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:19:02
	% EndTime: 2020-11-04 21:19:02
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t107 = cos(pkin(11));
	t106 = sin(pkin(11));
	t1 = [t107, -t106, 0, 0; t106, t107, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:19:02
	% EndTime: 2020-11-04 21:19:02
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->19), div. (0->0), fcn. (36->6), ass. (0->11)
	t108 = sin(pkin(11));
	t109 = sin(pkin(6));
	t117 = t108 * t109;
	t110 = cos(pkin(11));
	t116 = t110 * t109;
	t111 = cos(pkin(6));
	t112 = sin(qJ(2));
	t115 = t111 * t112;
	t113 = cos(qJ(2));
	t114 = t111 * t113;
	t1 = [-t108 * t115 + t110 * t113, -t108 * t114 - t110 * t112, t117, t110 * pkin(1) + pkin(7) * t117 + 0; t108 * t113 + t110 * t115, -t108 * t112 + t110 * t114, -t116, t108 * pkin(1) - pkin(7) * t116 + 0; t109 * t112, t109 * t113, t111, t111 * pkin(7) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:19:02
	% EndTime: 2020-11-04 21:19:02
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (29->27), mult. (64->48), div. (0->0), fcn. (85->8), ass. (0->18)
	t121 = sin(pkin(6));
	t134 = pkin(7) * t121;
	t120 = sin(pkin(11));
	t123 = cos(pkin(6));
	t133 = t120 * t123;
	t124 = sin(qJ(3));
	t132 = t121 * t124;
	t126 = cos(qJ(3));
	t131 = t121 * t126;
	t122 = cos(pkin(11));
	t130 = t122 * t123;
	t125 = sin(qJ(2));
	t129 = t123 * t125;
	t127 = cos(qJ(2));
	t128 = t123 * t127;
	t119 = t120 * t127 + t122 * t129;
	t118 = t120 * t129 - t122 * t127;
	t1 = [-t118 * t126 + t120 * t132, t118 * t124 + t120 * t131, t120 * t128 + t122 * t125, (pkin(2) * t122 + pkin(8) * t133) * t127 + (-pkin(2) * t133 + pkin(8) * t122) * t125 + t120 * t134 + t122 * pkin(1) + 0; t119 * t126 - t122 * t132, -t119 * t124 - t122 * t131, t120 * t125 - t122 * t128, (pkin(2) * t120 - pkin(8) * t130) * t127 + (pkin(2) * t130 + pkin(8) * t120) * t125 - t122 * t134 + t120 * pkin(1) + 0; t123 * t124 + t125 * t131, t123 * t126 - t125 * t132, -t121 * t127, t123 * pkin(7) + qJ(1) + 0 + (pkin(2) * t125 - pkin(8) * t127) * t121; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:19:02
	% EndTime: 2020-11-04 21:19:02
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (57->32), mult. (136->56), div. (0->0), fcn. (170->10), ass. (0->29)
	t142 = sin(pkin(6));
	t144 = cos(pkin(6));
	t147 = sin(qJ(2));
	t150 = cos(qJ(2));
	t146 = sin(qJ(3));
	t149 = cos(qJ(3));
	t154 = t149 * pkin(3) + t146 * pkin(9) + pkin(2);
	t152 = -pkin(8) * t150 + t154 * t147;
	t153 = t146 * pkin(3) - t149 * pkin(9) + pkin(7);
	t165 = t153 * t142 - t152 * t144;
	t161 = t142 * t149;
	t160 = t142 * t150;
	t159 = t144 * t147;
	t158 = t144 * t150;
	t157 = t146 * t142;
	t156 = t147 * t149;
	t155 = t149 * t150;
	t151 = pkin(8) * t147 + t154 * t150 + pkin(1);
	t148 = cos(qJ(4));
	t145 = sin(qJ(4));
	t143 = cos(pkin(11));
	t141 = sin(pkin(11));
	t140 = t144 * t156 - t157;
	t139 = t142 * t156 + t144 * t146;
	t138 = t141 * t158 + t143 * t147;
	t137 = t141 * t147 - t143 * t158;
	t136 = -t141 * t140 + t143 * t155;
	t135 = t143 * t140 + t141 * t155;
	t1 = [t136 * t148 + t138 * t145, -t136 * t145 + t138 * t148, -(t141 * t159 - t143 * t150) * t146 - t141 * t161, t165 * t141 + t151 * t143 + 0; t135 * t148 + t137 * t145, -t135 * t145 + t137 * t148, (t141 * t150 + t143 * t159) * t146 + t143 * t161, t151 * t141 - t165 * t143 + 0; t139 * t148 - t145 * t160, -t139 * t145 - t148 * t160, -t144 * t149 + t147 * t157, t152 * t142 + t153 * t144 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:19:02
	% EndTime: 2020-11-04 21:19:02
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (80->45), mult. (150->70), div. (0->0), fcn. (196->12), ass. (0->33)
	t198 = pkin(4) * sin(qJ(4));
	t181 = sin(pkin(6));
	t197 = t181 * pkin(7);
	t180 = sin(pkin(11));
	t183 = cos(pkin(6));
	t196 = t180 * t183;
	t185 = sin(qJ(3));
	t195 = t181 * t185;
	t187 = cos(qJ(3));
	t194 = t181 * t187;
	t188 = cos(qJ(2));
	t193 = t181 * t188;
	t182 = cos(pkin(11));
	t192 = t182 * t183;
	t186 = sin(qJ(2));
	t191 = t183 * t186;
	t190 = t183 * t188;
	t189 = -pkin(10) - pkin(9);
	t179 = qJ(4) + qJ(5);
	t178 = cos(t179);
	t177 = sin(t179);
	t176 = cos(qJ(4)) * pkin(4) + pkin(3);
	t175 = t183 * t185 + t186 * t194;
	t174 = -t183 * t187 + t186 * t195;
	t173 = t180 * t190 + t182 * t186;
	t172 = t180 * t188 + t182 * t191;
	t171 = t180 * t186 - t182 * t190;
	t170 = t180 * t191 - t182 * t188;
	t169 = -t170 * t187 + t180 * t195;
	t168 = t172 * t187 - t182 * t195;
	t167 = t172 * t185 + t182 * t194;
	t166 = t170 * t185 + t180 * t194;
	t1 = [t169 * t178 + t173 * t177, -t169 * t177 + t173 * t178, -t166, t169 * t176 + t166 * t189 + t173 * t198 + (t182 * pkin(2) + pkin(8) * t196) * t188 + (-pkin(2) * t196 + t182 * pkin(8)) * t186 + t180 * t197 + t182 * pkin(1) + 0; t168 * t178 + t171 * t177, -t168 * t177 + t171 * t178, t167, t168 * t176 - t167 * t189 + t171 * t198 + (t180 * pkin(2) - pkin(8) * t192) * t188 + (pkin(2) * t192 + t180 * pkin(8)) * t186 - t182 * t197 + t180 * pkin(1) + 0; t175 * t178 - t177 * t193, -t175 * t177 - t178 * t193, t174, t183 * pkin(7) - t174 * t189 + t175 * t176 + qJ(1) + 0 + (pkin(2) * t186 + (-pkin(8) - t198) * t188) * t181; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:19:02
	% EndTime: 2020-11-04 21:19:02
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (118->51), mult. (212->76), div. (0->0), fcn. (278->12), ass. (0->39)
	t237 = pkin(4) * sin(qJ(4));
	t220 = sin(pkin(6));
	t236 = t220 * pkin(7);
	t219 = sin(pkin(11));
	t222 = cos(pkin(6));
	t235 = t219 * t222;
	t224 = sin(qJ(3));
	t234 = t220 * t224;
	t226 = cos(qJ(3));
	t233 = t220 * t226;
	t227 = cos(qJ(2));
	t232 = t220 * t227;
	t221 = cos(pkin(11));
	t231 = t221 * t222;
	t225 = sin(qJ(2));
	t230 = t222 * t225;
	t229 = t222 * t227;
	t228 = -pkin(10) - pkin(9);
	t218 = qJ(4) + qJ(5);
	t217 = cos(t218);
	t216 = sin(t218);
	t215 = cos(qJ(4)) * pkin(4) + pkin(3);
	t214 = t222 * t224 + t225 * t233;
	t213 = -t222 * t226 + t225 * t234;
	t212 = t219 * t229 + t221 * t225;
	t211 = t219 * t227 + t221 * t230;
	t210 = t219 * t225 - t221 * t229;
	t209 = t219 * t230 - t221 * t227;
	t208 = -t209 * t226 + t219 * t234;
	t207 = t211 * t226 - t221 * t234;
	t206 = t211 * t224 + t221 * t233;
	t205 = t209 * t224 + t219 * t233;
	t204 = t214 * t217 - t216 * t232;
	t203 = t214 * t216 + t217 * t232;
	t202 = t208 * t217 + t212 * t216;
	t201 = t208 * t216 - t212 * t217;
	t200 = t207 * t217 + t210 * t216;
	t199 = t207 * t216 - t210 * t217;
	t1 = [t202, -t205, t201, t202 * pkin(5) + t201 * qJ(6) + t208 * t215 + t205 * t228 + t212 * t237 + (t221 * pkin(2) + pkin(8) * t235) * t227 + (-pkin(2) * t235 + t221 * pkin(8)) * t225 + t219 * t236 + t221 * pkin(1) + 0; t200, t206, t199, t200 * pkin(5) + t199 * qJ(6) + t207 * t215 - t206 * t228 + t210 * t237 + (t219 * pkin(2) - pkin(8) * t231) * t227 + (pkin(2) * t231 + t219 * pkin(8)) * t225 - t221 * t236 + t219 * pkin(1) + 0; t204, t213, t203, t204 * pkin(5) + t222 * pkin(7) + t203 * qJ(6) - t213 * t228 + t214 * t215 + qJ(1) + 0 + (pkin(2) * t225 + (-pkin(8) - t237) * t227) * t220; 0, 0, 0, 1;];
	Tc_mdh = t1;
end