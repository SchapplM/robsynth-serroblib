% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRPPR9 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:25
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRRPPR9_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR9_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPPR9_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR9_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:25:14
	% EndTime: 2020-11-04 22:25:14
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:25:14
	% EndTime: 2020-11-04 22:25:14
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t108 = cos(qJ(1));
	t107 = sin(qJ(1));
	t1 = [t108, -t107, 0, 0; t107, t108, 0, 0; 0, 0, 1, pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:25:14
	% EndTime: 2020-11-04 22:25:14
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->17), div. (0->0), fcn. (36->6), ass. (0->13)
	t109 = sin(pkin(6));
	t112 = sin(qJ(1));
	t120 = t112 * t109;
	t111 = sin(qJ(2));
	t119 = t112 * t111;
	t113 = cos(qJ(2));
	t118 = t112 * t113;
	t114 = cos(qJ(1));
	t117 = t114 * t109;
	t116 = t114 * t111;
	t115 = t114 * t113;
	t110 = cos(pkin(6));
	t1 = [-t110 * t119 + t115, -t110 * t118 - t116, t120, t114 * pkin(1) + pkin(8) * t120 + 0; t110 * t116 + t118, t110 * t115 - t119, -t117, t112 * pkin(1) - pkin(8) * t117 + 0; t109 * t111, t109 * t113, t110, t110 * pkin(8) + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:25:14
	% EndTime: 2020-11-04 22:25:14
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (29->23), mult. (60->37), div. (0->0), fcn. (81->8), ass. (0->20)
	t124 = sin(pkin(6));
	t129 = cos(qJ(3));
	t139 = t124 * t129;
	t126 = sin(qJ(3));
	t138 = t126 * t124;
	t127 = sin(qJ(2));
	t137 = t127 * t129;
	t128 = sin(qJ(1));
	t136 = t128 * t127;
	t130 = cos(qJ(2));
	t135 = t128 * t130;
	t131 = cos(qJ(1));
	t134 = t131 * t127;
	t133 = t131 * t130;
	t132 = pkin(2) * t127 - pkin(9) * t130;
	t125 = cos(pkin(6));
	t123 = t130 * pkin(2) + t127 * pkin(9) + pkin(1);
	t122 = t125 * t134 + t135;
	t121 = t124 * pkin(8) - t132 * t125;
	t1 = [(-t125 * t137 + t138) * t128 + t129 * t133, (t125 * t136 - t133) * t126 + t128 * t139, t125 * t135 + t134, t121 * t128 + t123 * t131 + 0; t122 * t129 - t131 * t138, -t122 * t126 - t131 * t139, -t125 * t133 + t136, -t121 * t131 + t123 * t128 + 0; t124 * t137 + t125 * t126, t125 * t129 - t127 * t138, -t124 * t130, t125 * pkin(8) + t132 * t124 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:25:14
	% EndTime: 2020-11-04 22:25:14
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (57->32), mult. (112->56), div. (0->0), fcn. (146->10), ass. (0->29)
	t150 = sin(pkin(6));
	t156 = cos(qJ(3));
	t167 = t150 * t156;
	t157 = cos(qJ(2));
	t166 = t150 * t157;
	t152 = cos(pkin(6));
	t154 = sin(qJ(2));
	t165 = t152 * t154;
	t164 = t152 * t157;
	t153 = sin(qJ(3));
	t163 = t153 * t150;
	t162 = t154 * t156;
	t161 = t156 * t157;
	t148 = t156 * pkin(3) + qJ(4) * t153 + pkin(2);
	t160 = t157 * pkin(9) - t148 * t154;
	t159 = t153 * pkin(3) - qJ(4) * t156 + pkin(8);
	t158 = cos(qJ(1));
	t155 = sin(qJ(1));
	t151 = cos(pkin(11));
	t149 = sin(pkin(11));
	t147 = t149 * t154 + t151 * t161;
	t146 = t152 * t162 - t163;
	t145 = t150 * t162 + t152 * t153;
	t144 = t149 * t161 - t154 * t151;
	t143 = t154 * pkin(9) + t148 * t157 + pkin(1);
	t142 = t149 * t146 + t151 * t164;
	t141 = -t151 * t146 + t149 * t164;
	t140 = t150 * t159 + t160 * t152;
	t1 = [t141 * t155 + t158 * t147, t142 * t155 - t158 * t144, (-t155 * t165 + t158 * t157) * t153 - t155 * t167, t140 * t155 + t143 * t158 + 0; -t141 * t158 + t155 * t147, -t142 * t158 - t155 * t144, (t155 * t157 + t158 * t165) * t153 + t158 * t167, -t140 * t158 + t143 * t155 + 0; t145 * t151 - t149 * t166, -t145 * t149 - t151 * t166, -t152 * t156 + t154 * t163, -t160 * t150 + t159 * t152 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:25:14
	% EndTime: 2020-11-04 22:25:14
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (83->36), mult. (138->60), div. (0->0), fcn. (172->10), ass. (0->31)
	t180 = sin(pkin(6));
	t186 = cos(qJ(3));
	t197 = t180 * t186;
	t187 = cos(qJ(2));
	t196 = t180 * t187;
	t182 = cos(pkin(6));
	t184 = sin(qJ(2));
	t195 = t182 * t184;
	t194 = t182 * t187;
	t183 = sin(qJ(3));
	t193 = t183 * t180;
	t192 = t184 * t186;
	t191 = t186 * t187;
	t179 = sin(pkin(11));
	t181 = cos(pkin(11));
	t177 = pkin(4) * t181 + qJ(5) * t179 + pkin(3);
	t172 = qJ(4) * t183 + t177 * t186 + pkin(2);
	t178 = -t179 * pkin(4) + qJ(5) * t181 - pkin(9);
	t190 = t172 * t184 + t178 * t187;
	t189 = qJ(4) * t186 - t177 * t183 - pkin(8);
	t188 = cos(qJ(1));
	t185 = sin(qJ(1));
	t176 = t179 * t184 + t181 * t191;
	t175 = t182 * t192 - t193;
	t174 = t180 * t192 + t182 * t183;
	t173 = t179 * t191 - t184 * t181;
	t171 = t179 * t175 + t181 * t194;
	t170 = -t181 * t175 + t179 * t194;
	t169 = t172 * t187 - t178 * t184 + pkin(1);
	t168 = t189 * t180 + t190 * t182;
	t1 = [t170 * t185 + t188 * t176, (-t185 * t195 + t188 * t187) * t183 - t185 * t197, -t171 * t185 + t188 * t173, -t168 * t185 + t169 * t188 + 0; -t170 * t188 + t185 * t176, (t185 * t187 + t188 * t195) * t183 + t188 * t197, t171 * t188 + t185 * t173, t168 * t188 + t169 * t185 + 0; t174 * t181 - t179 * t196, -t182 * t186 + t184 * t193, t174 * t179 + t181 * t196, t190 * t180 - t189 * t182 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:25:14
	% EndTime: 2020-11-04 22:25:15
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (130->44), mult. (206->72), div. (0->0), fcn. (266->12), ass. (0->41)
	t213 = sin(pkin(11));
	t215 = cos(pkin(11));
	t226 = pkin(4) + pkin(5);
	t211 = qJ(5) * t213 + t226 * t215 + pkin(3);
	t217 = qJ(4) - pkin(10);
	t219 = sin(qJ(3));
	t223 = cos(qJ(3));
	t240 = -t211 * t219 + t217 * t223 - pkin(8);
	t214 = sin(pkin(6));
	t238 = t214 * t223;
	t224 = cos(qJ(2));
	t237 = t214 * t224;
	t216 = cos(pkin(6));
	t220 = sin(qJ(2));
	t236 = t216 * t220;
	t235 = t216 * t224;
	t234 = t219 * t214;
	t233 = t220 * t223;
	t232 = t223 * t224;
	t227 = -t216 * t233 + t234;
	t202 = t213 * t235 + t215 * t227;
	t203 = -t213 * t227 + t215 * t235;
	t218 = sin(qJ(6));
	t222 = cos(qJ(6));
	t230 = t202 * t222 - t218 * t203;
	t229 = t202 * t218 + t203 * t222;
	t206 = t211 * t223 + t217 * t219 + pkin(2);
	t212 = qJ(5) * t215 - t226 * t213 - pkin(9);
	t228 = t206 * t220 + t212 * t224;
	t225 = cos(qJ(1));
	t221 = sin(qJ(1));
	t210 = t213 * t220 + t215 * t232;
	t209 = t214 * t233 + t216 * t219;
	t208 = t213 * t232 - t220 * t215;
	t205 = t209 * t215 - t213 * t237;
	t204 = t209 * t213 + t215 * t237;
	t201 = t218 * t208 + t210 * t222;
	t200 = t208 * t222 - t218 * t210;
	t199 = t206 * t224 - t212 * t220 + pkin(1);
	t198 = t214 * t240 + t228 * t216;
	t1 = [t201 * t225 + t230 * t221, t200 * t225 - t229 * t221, (t221 * t236 - t225 * t224) * t219 + t221 * t238, -t198 * t221 + t199 * t225 + 0; t201 * t221 - t230 * t225, t200 * t221 + t229 * t225, (-t221 * t224 - t225 * t236) * t219 - t225 * t238, t198 * t225 + t199 * t221 + 0; t204 * t218 + t205 * t222, t204 * t222 - t218 * t205, t216 * t223 - t220 * t234, t228 * t214 - t240 * t216 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end