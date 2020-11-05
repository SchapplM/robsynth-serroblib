% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPRPR6 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:09
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRPRPR6_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR6_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPR6_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR6_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:09:48
	% EndTime: 2020-11-04 22:09:48
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:09:48
	% EndTime: 2020-11-04 22:09:48
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t112 = cos(qJ(1));
	t111 = sin(qJ(1));
	t1 = [t112, -t111, 0, 0; t111, t112, 0, 0; 0, 0, 1, pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:09:48
	% EndTime: 2020-11-04 22:09:49
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (11->11), mult. (23->17), div. (0->0), fcn. (36->6), ass. (0->13)
	t115 = sin(qJ(2));
	t118 = cos(qJ(1));
	t124 = t115 * t118;
	t113 = sin(pkin(6));
	t116 = sin(qJ(1));
	t123 = t116 * t113;
	t122 = t116 * t115;
	t117 = cos(qJ(2));
	t121 = t116 * t117;
	t120 = t117 * t118;
	t119 = t118 * t113;
	t114 = cos(pkin(6));
	t1 = [-t114 * t122 + t120, -t114 * t121 - t124, t123, pkin(1) * t118 + pkin(8) * t123 + 0; t114 * t124 + t121, t114 * t120 - t122, -t119, t116 * pkin(1) - pkin(8) * t119 + 0; t113 * t115, t113 * t117, t114, pkin(8) * t114 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:09:49
	% EndTime: 2020-11-04 22:09:49
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (53->22), mult. (33->28), div. (0->0), fcn. (44->12), ass. (0->20)
	t144 = pkin(2) * sin(qJ(2));
	t137 = qJ(2) + pkin(11);
	t143 = cos(qJ(1));
	t142 = sin(qJ(1));
	t140 = pkin(8) + qJ(3);
	t139 = cos(pkin(6));
	t138 = sin(pkin(6));
	t136 = cos(t137);
	t135 = sin(t137);
	t134 = pkin(6) - t137;
	t133 = pkin(6) + t137;
	t132 = cos(qJ(2)) * pkin(2) + pkin(1);
	t131 = cos(t134);
	t130 = cos(t133);
	t129 = sin(t134);
	t128 = sin(t133);
	t127 = t130 + t131;
	t126 = -t128 + t129;
	t125 = -t138 * t140 + t139 * t144;
	t1 = [t143 * t136 + t142 * t126 / 0.2e1, -t143 * t135 - t142 * t127 / 0.2e1, t142 * t138, -t125 * t142 + t143 * t132 + 0; t142 * t136 - t143 * t126 / 0.2e1, -t142 * t135 + t143 * t127 / 0.2e1, -t143 * t138, t143 * t125 + t142 * t132 + 0; t131 / 0.2e1 - t130 / 0.2e1, t129 / 0.2e1 + t128 / 0.2e1, t139, t138 * t144 + t139 * t140 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:09:49
	% EndTime: 2020-11-04 22:09:49
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (73->30), mult. (81->43), div. (0->0), fcn. (101->16), ass. (0->30)
	t158 = sin(pkin(6));
	t162 = sin(qJ(4));
	t173 = t158 * t162;
	t165 = cos(qJ(4));
	t172 = t158 * t165;
	t167 = cos(qJ(1));
	t171 = t158 * t167;
	t156 = qJ(2) + pkin(11);
	t154 = sin(t156);
	t164 = sin(qJ(1));
	t170 = t164 * t154;
	t169 = t167 * t154;
	t157 = sin(pkin(11));
	t159 = cos(pkin(11));
	t150 = t159 * pkin(3) + t157 * pkin(9) + pkin(2);
	t151 = -t157 * pkin(3) + t159 * pkin(9);
	t163 = sin(qJ(2));
	t166 = cos(qJ(2));
	t168 = t150 * t163 - t151 * t166;
	t161 = pkin(8) + qJ(3);
	t160 = cos(pkin(6));
	t155 = cos(t156);
	t153 = pkin(6) - t156;
	t152 = pkin(6) + t156;
	t149 = cos(t152) + cos(t153);
	t148 = t164 * t155 + t160 * t169;
	t147 = -t167 * t155 + t160 * t170;
	t146 = t150 * t166 + t151 * t163 + pkin(1);
	t145 = t158 * t161 - t168 * t160;
	t1 = [-t147 * t165 + t164 * t173, t147 * t162 + t164 * t172, t169 + t164 * t149 / 0.2e1, t145 * t164 + t146 * t167 + 0; t148 * t165 - t162 * t171, -t148 * t162 - t165 * t171, t170 - t167 * t149 / 0.2e1, -t145 * t167 + t146 * t164 + 0; t154 * t172 + t160 * t162, -t154 * t173 + t160 * t165, -sin(t153) / 0.2e1 - sin(t152) / 0.2e1, t168 * t158 + t160 * t161 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:09:49
	% EndTime: 2020-11-04 22:09:49
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (99->40), mult. (117->49), div. (0->0), fcn. (147->16), ass. (0->36)
	t193 = sin(pkin(6));
	t197 = sin(qJ(4));
	t208 = t193 * t197;
	t200 = cos(qJ(4));
	t207 = t193 * t200;
	t202 = cos(qJ(1));
	t206 = t193 * t202;
	t191 = qJ(2) + pkin(11);
	t189 = sin(t191);
	t199 = sin(qJ(1));
	t205 = t199 * t189;
	t204 = t202 * t189;
	t192 = sin(pkin(11));
	t194 = cos(pkin(11));
	t185 = t194 * pkin(3) + t192 * pkin(9) + pkin(2);
	t186 = -t192 * pkin(3) + t194 * pkin(9);
	t198 = sin(qJ(2));
	t201 = cos(qJ(2));
	t203 = t185 * t198 - t186 * t201;
	t196 = pkin(8) + qJ(3);
	t195 = cos(pkin(6));
	t190 = cos(t191);
	t188 = pkin(6) - t191;
	t187 = pkin(6) + t191;
	t184 = cos(t187) + cos(t188);
	t183 = t189 * t207 + t195 * t197;
	t182 = t189 * t208 - t195 * t200;
	t181 = t199 * t190 + t195 * t204;
	t180 = -t202 * t190 + t195 * t205;
	t179 = t185 * t201 + t186 * t198 + pkin(1);
	t178 = -t180 * t200 + t199 * t208;
	t177 = t181 * t200 - t197 * t206;
	t176 = t181 * t197 + t200 * t206;
	t175 = t180 * t197 + t199 * t207;
	t174 = t193 * t196 - t203 * t195;
	t1 = [t204 + t199 * t184 / 0.2e1, -t178, -t175, t178 * pkin(4) - t175 * qJ(5) + t174 * t199 + t179 * t202 + 0; t205 - t202 * t184 / 0.2e1, -t177, t176, t177 * pkin(4) + t176 * qJ(5) - t174 * t202 + t179 * t199 + 0; -sin(t188) / 0.2e1 - sin(t187) / 0.2e1, -t183, t182, t183 * pkin(4) + t182 * qJ(5) + t203 * t193 + t195 * t196 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:09:49
	% EndTime: 2020-11-04 22:09:49
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (127->55), mult. (174->80), div. (0->0), fcn. (203->16), ass. (0->41)
	t217 = sin(pkin(11));
	t219 = cos(pkin(11));
	t229 = pkin(5) + pkin(9);
	t222 = sin(qJ(4));
	t226 = cos(qJ(4));
	t230 = pkin(4) + pkin(10);
	t233 = qJ(5) * t222 + t230 * t226;
	t251 = t217 * (pkin(3) + t233) - t229 * t219;
	t250 = t219 * pkin(3) + pkin(2);
	t216 = qJ(2) + pkin(11);
	t214 = cos(t216);
	t218 = sin(pkin(6));
	t249 = t214 * t218;
	t248 = t218 * t222;
	t220 = cos(pkin(6));
	t247 = t220 * t226;
	t221 = sin(qJ(6));
	t224 = sin(qJ(1));
	t246 = t224 * t221;
	t225 = cos(qJ(6));
	t245 = t225 * t224;
	t244 = t226 * t218;
	t228 = cos(qJ(1));
	t243 = t226 * t228;
	t242 = t228 * t221;
	t241 = t228 * t225;
	t239 = t222 * t246;
	t238 = t222 * t245;
	t237 = t224 * t244;
	t236 = t218 * t243;
	t235 = t222 * t242;
	t234 = t222 * t241;
	t231 = -qJ(5) * t226 + t230 * t222 + pkin(8) + qJ(3);
	t227 = cos(qJ(2));
	t223 = sin(qJ(2));
	t213 = sin(t216);
	t212 = t213 * t248 - t247;
	t211 = t229 * t217 + t233 * t219 + t250;
	t210 = t211 * t227 - t251 * t223 + pkin(1);
	t209 = -t231 * t218 + (t211 * t223 + t251 * t227) * t220;
	t1 = [(t220 * t245 + t235) * t214 + (-t220 * t239 + t241) * t213 - t221 * t237, (-t220 * t246 + t234) * t214 + (-t220 * t238 - t242) * t213 - t225 * t237, t214 * t243 + (-t213 * t247 + t248) * t224, -t209 * t224 + t210 * t228 + 0; (-t220 * t241 + t239) * t214 + (t220 * t235 + t245) * t213 + t221 * t236, (t220 * t242 + t238) * t214 + (t220 * t234 - t246) * t213 + t225 * t236, (t228 * t220 * t213 + t224 * t214) * t226 - t228 * t248, t209 * t228 + t210 * t224 + 0; t212 * t221 - t225 * t249, t212 * t225 + t221 * t249, t213 * t244 + t220 * t222, pkin(7) + 0 + (-sin(pkin(6) - t216) / 0.2e1 - sin(pkin(6) + t216) / 0.2e1) * pkin(5) + t231 * t220 + (t233 * t213 + (t217 * pkin(9) + t250) * t223 - (-t217 * pkin(3) + t219 * pkin(9)) * t227) * t218; 0, 0, 0, 1;];
	Tc_mdh = t1;
end