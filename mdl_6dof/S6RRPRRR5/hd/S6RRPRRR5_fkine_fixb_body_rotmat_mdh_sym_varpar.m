% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPRRR5 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:18
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRPRRR5_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR5_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRR5_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR5_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [12x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:18:36
	% EndTime: 2020-11-04 22:18:36
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:18:36
	% EndTime: 2020-11-04 22:18:36
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t116 = cos(qJ(1));
	t115 = sin(qJ(1));
	t1 = [t116, -t115, 0, 0; t115, t116, 0, 0; 0, 0, 1, pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:18:36
	% EndTime: 2020-11-04 22:18:36
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->17), div. (0->0), fcn. (36->6), ass. (0->13)
	t117 = sin(pkin(6));
	t120 = sin(qJ(1));
	t128 = t120 * t117;
	t119 = sin(qJ(2));
	t127 = t120 * t119;
	t121 = cos(qJ(2));
	t126 = t120 * t121;
	t122 = cos(qJ(1));
	t125 = t122 * t117;
	t124 = t122 * t119;
	t123 = t122 * t121;
	t118 = cos(pkin(6));
	t1 = [-t118 * t127 + t123, -t118 * t126 - t124, t128, t122 * pkin(1) + pkin(8) * t128 + 0; t118 * t124 + t126, t118 * t123 - t127, -t125, t120 * pkin(1) - pkin(8) * t125 + 0; t117 * t119, t117 * t121, t118, t118 * pkin(8) + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:18:36
	% EndTime: 2020-11-04 22:18:36
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (53->22), mult. (33->28), div. (0->0), fcn. (44->12), ass. (0->20)
	t148 = pkin(2) * sin(qJ(2));
	t141 = qJ(2) + pkin(12);
	t147 = cos(qJ(1));
	t146 = sin(qJ(1));
	t144 = pkin(8) + qJ(3);
	t143 = cos(pkin(6));
	t142 = sin(pkin(6));
	t140 = cos(t141);
	t139 = sin(t141);
	t138 = pkin(6) - t141;
	t137 = pkin(6) + t141;
	t136 = cos(qJ(2)) * pkin(2) + pkin(1);
	t135 = cos(t138);
	t134 = cos(t137);
	t133 = sin(t138);
	t132 = sin(t137);
	t131 = t134 + t135;
	t130 = -t132 + t133;
	t129 = -t142 * t144 + t143 * t148;
	t1 = [t147 * t140 + t146 * t130 / 0.2e1, -t147 * t139 - t146 * t131 / 0.2e1, t146 * t142, -t129 * t146 + t147 * t136 + 0; t146 * t140 - t147 * t130 / 0.2e1, -t146 * t139 + t147 * t131 / 0.2e1, -t147 * t142, t147 * t129 + t146 * t136 + 0; t135 / 0.2e1 - t134 / 0.2e1, t133 / 0.2e1 + t132 / 0.2e1, t143, t142 * t148 + t143 * t144 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:18:36
	% EndTime: 2020-11-04 22:18:36
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (73->30), mult. (81->43), div. (0->0), fcn. (101->16), ass. (0->30)
	t162 = sin(pkin(6));
	t166 = sin(qJ(4));
	t177 = t162 * t166;
	t169 = cos(qJ(4));
	t176 = t162 * t169;
	t171 = cos(qJ(1));
	t175 = t162 * t171;
	t160 = qJ(2) + pkin(12);
	t158 = sin(t160);
	t168 = sin(qJ(1));
	t174 = t168 * t158;
	t173 = t171 * t158;
	t161 = sin(pkin(12));
	t163 = cos(pkin(12));
	t154 = t163 * pkin(3) + t161 * pkin(9) + pkin(2);
	t155 = -t161 * pkin(3) + t163 * pkin(9);
	t167 = sin(qJ(2));
	t170 = cos(qJ(2));
	t172 = t154 * t167 - t155 * t170;
	t165 = pkin(8) + qJ(3);
	t164 = cos(pkin(6));
	t159 = cos(t160);
	t157 = pkin(6) - t160;
	t156 = pkin(6) + t160;
	t153 = cos(t156) + cos(t157);
	t152 = t168 * t159 + t164 * t173;
	t151 = -t171 * t159 + t164 * t174;
	t150 = t154 * t170 + t155 * t167 + pkin(1);
	t149 = t162 * t165 - t172 * t164;
	t1 = [-t151 * t169 + t168 * t177, t151 * t166 + t168 * t176, t173 + t168 * t153 / 0.2e1, t149 * t168 + t150 * t171 + 0; t152 * t169 - t166 * t175, -t152 * t166 - t169 * t175, t174 - t171 * t153 / 0.2e1, -t149 * t171 + t150 * t168 + 0; t158 * t176 + t164 * t166, -t158 * t177 + t164 * t169, -sin(t157) / 0.2e1 - sin(t156) / 0.2e1, t172 * t162 + t164 * t165 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:18:36
	% EndTime: 2020-11-04 22:18:36
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (102->44), mult. (170->73), div. (0->0), fcn. (204->14), ass. (0->39)
	t185 = qJ(2) + pkin(12);
	t183 = sin(t185);
	t189 = cos(pkin(6));
	t215 = t183 * t189;
	t187 = sin(pkin(6));
	t191 = sin(qJ(4));
	t214 = t187 * t191;
	t194 = cos(qJ(5));
	t213 = t187 * t194;
	t190 = sin(qJ(5));
	t193 = sin(qJ(1));
	t212 = t190 * t193;
	t197 = cos(qJ(1));
	t211 = t190 * t197;
	t195 = cos(qJ(4));
	t210 = t193 * t195;
	t209 = t194 * t193;
	t208 = t194 * t197;
	t207 = t195 * t197;
	t206 = t190 * t214;
	t205 = t191 * t213;
	t204 = t195 * t209;
	t203 = t190 * t210;
	t202 = t190 * t207;
	t201 = t194 * t207;
	t186 = sin(pkin(12));
	t188 = cos(pkin(12));
	t199 = pkin(4) * t195 + pkin(10) * t191 + pkin(3);
	t180 = t186 * pkin(9) + t199 * t188 + pkin(2);
	t181 = t188 * pkin(9) - t199 * t186;
	t192 = sin(qJ(2));
	t196 = cos(qJ(2));
	t200 = t180 * t192 - t181 * t196;
	t198 = t191 * pkin(4) - pkin(10) * t195 + pkin(8) + qJ(3);
	t184 = cos(t185);
	t182 = t187 * t183 * t195 + t189 * t191;
	t179 = t180 * t196 + t181 * t192 + pkin(1);
	t178 = t187 * t198 - t200 * t189;
	t1 = [(t189 * t212 + t201) * t184 + (-t189 * t204 + t211) * t183 + t193 * t205, (t189 * t209 - t202) * t184 + (t189 * t203 + t208) * t183 - t193 * t206, -(-t197 * t184 + t193 * t215) * t191 - t187 * t210, t178 * t193 + t179 * t197 + 0; (-t189 * t211 + t204) * t184 + (t189 * t201 + t212) * t183 - t197 * t205, (-t189 * t208 - t203) * t184 + (-t189 * t202 + t209) * t183 + t197 * t206, (t193 * t184 + t197 * t215) * t191 + t187 * t207, -t178 * t197 + t179 * t193 + 0; -t187 * t184 * t190 + t182 * t194, -t182 * t190 - t184 * t213, t183 * t214 - t189 * t195, t200 * t187 + t198 * t189 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:18:36
	% EndTime: 2020-11-04 22:18:36
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (174->49), mult. (169->66), div. (0->0), fcn. (212->20), ass. (0->45)
	t260 = pkin(5) * sin(qJ(5));
	t242 = sin(pkin(6));
	t247 = sin(qJ(4));
	t259 = t242 * t247;
	t250 = cos(qJ(4));
	t258 = t242 * t250;
	t252 = cos(qJ(1));
	t257 = t242 * t252;
	t239 = qJ(2) + pkin(12);
	t235 = sin(t239);
	t249 = sin(qJ(1));
	t256 = t249 * t235;
	t255 = t252 * t235;
	t241 = sin(pkin(12));
	t243 = cos(pkin(12));
	t230 = t243 * pkin(3) + t241 * pkin(9) + pkin(2);
	t231 = -t241 * pkin(3) + t243 * pkin(9);
	t248 = sin(qJ(2));
	t251 = cos(qJ(2));
	t254 = t230 * t248 - t231 * t251;
	t253 = -pkin(11) - pkin(10);
	t245 = pkin(8) + qJ(3);
	t244 = cos(pkin(6));
	t240 = qJ(5) + qJ(6);
	t238 = cos(t240);
	t237 = sin(t240);
	t236 = cos(t239);
	t234 = pkin(6) - t239;
	t233 = pkin(6) + t239;
	t232 = cos(qJ(5)) * pkin(5) + pkin(4);
	t229 = cos(t233) + cos(t234);
	t228 = -sin(t234) / 0.2e1 - sin(t233) / 0.2e1;
	t227 = t235 * t258 + t244 * t247;
	t226 = t235 * t259 - t244 * t250;
	t225 = t249 * t236 + t244 * t255;
	t224 = -t252 * t236 + t244 * t256;
	t223 = t255 + t249 * t229 / 0.2e1;
	t222 = t256 - t252 * t229 / 0.2e1;
	t221 = t230 * t251 + t231 * t248 + pkin(1);
	t220 = -t224 * t250 + t249 * t259;
	t219 = t225 * t250 - t247 * t257;
	t218 = t225 * t247 + t250 * t257;
	t217 = t224 * t247 + t249 * t258;
	t216 = t242 * t245 - t254 * t244;
	t1 = [t220 * t238 + t223 * t237, -t220 * t237 + t223 * t238, -t217, t216 * t249 + t217 * t253 + t220 * t232 + t221 * t252 + t223 * t260 + 0; t219 * t238 + t222 * t237, -t219 * t237 + t222 * t238, t218, -t216 * t252 - t218 * t253 + t219 * t232 + t221 * t249 + t222 * t260 + 0; t227 * t238 + t228 * t237, -t227 * t237 + t228 * t238, t226, -t226 * t253 + t227 * t232 + t228 * t260 + t254 * t242 + t244 * t245 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end