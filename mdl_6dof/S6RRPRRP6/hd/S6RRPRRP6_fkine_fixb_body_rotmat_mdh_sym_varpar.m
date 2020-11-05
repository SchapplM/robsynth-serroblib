% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPRRP6 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:14
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRPRRP6_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP6_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRP6_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP6_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:14:33
	% EndTime: 2020-11-04 22:14:33
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:14:33
	% EndTime: 2020-11-04 22:14:33
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t122 = cos(qJ(1));
	t121 = sin(qJ(1));
	t1 = [t122, -t121, 0, 0; t121, t122, 0, 0; 0, 0, 1, pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:14:33
	% EndTime: 2020-11-04 22:14:33
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->17), div. (0->0), fcn. (36->6), ass. (0->13)
	t123 = sin(pkin(6));
	t126 = sin(qJ(1));
	t134 = t126 * t123;
	t125 = sin(qJ(2));
	t133 = t126 * t125;
	t127 = cos(qJ(2));
	t132 = t126 * t127;
	t128 = cos(qJ(1));
	t131 = t128 * t123;
	t130 = t128 * t125;
	t129 = t128 * t127;
	t124 = cos(pkin(6));
	t1 = [-t124 * t133 + t129, -t124 * t132 - t130, t134, t128 * pkin(1) + pkin(8) * t134 + 0; t124 * t130 + t132, t124 * t129 - t133, -t131, t126 * pkin(1) - pkin(8) * t131 + 0; t123 * t125, t123 * t127, t124, t124 * pkin(8) + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:14:33
	% EndTime: 2020-11-04 22:14:33
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (53->22), mult. (33->28), div. (0->0), fcn. (44->12), ass. (0->20)
	t154 = pkin(2) * sin(qJ(2));
	t147 = qJ(2) + pkin(11);
	t153 = cos(qJ(1));
	t152 = sin(qJ(1));
	t150 = pkin(8) + qJ(3);
	t149 = cos(pkin(6));
	t148 = sin(pkin(6));
	t146 = cos(t147);
	t145 = sin(t147);
	t144 = pkin(6) - t147;
	t143 = pkin(6) + t147;
	t142 = cos(qJ(2)) * pkin(2) + pkin(1);
	t141 = cos(t144);
	t140 = cos(t143);
	t139 = sin(t144);
	t138 = sin(t143);
	t137 = t140 + t141;
	t136 = -t138 + t139;
	t135 = -t148 * t150 + t149 * t154;
	t1 = [t153 * t146 + t152 * t136 / 0.2e1, -t153 * t145 - t152 * t137 / 0.2e1, t152 * t148, -t135 * t152 + t153 * t142 + 0; t152 * t146 - t153 * t136 / 0.2e1, -t152 * t145 + t153 * t137 / 0.2e1, -t153 * t148, t153 * t135 + t152 * t142 + 0; t141 / 0.2e1 - t140 / 0.2e1, t139 / 0.2e1 + t138 / 0.2e1, t149, t148 * t154 + t149 * t150 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:14:33
	% EndTime: 2020-11-04 22:14:34
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (73->30), mult. (81->43), div. (0->0), fcn. (101->16), ass. (0->30)
	t168 = sin(pkin(6));
	t172 = sin(qJ(4));
	t183 = t168 * t172;
	t175 = cos(qJ(4));
	t182 = t168 * t175;
	t177 = cos(qJ(1));
	t181 = t168 * t177;
	t166 = qJ(2) + pkin(11);
	t164 = sin(t166);
	t174 = sin(qJ(1));
	t180 = t174 * t164;
	t179 = t177 * t164;
	t167 = sin(pkin(11));
	t169 = cos(pkin(11));
	t160 = t169 * pkin(3) + t167 * pkin(9) + pkin(2);
	t161 = -t167 * pkin(3) + t169 * pkin(9);
	t173 = sin(qJ(2));
	t176 = cos(qJ(2));
	t178 = t160 * t173 - t161 * t176;
	t171 = pkin(8) + qJ(3);
	t170 = cos(pkin(6));
	t165 = cos(t166);
	t163 = pkin(6) - t166;
	t162 = pkin(6) + t166;
	t159 = cos(t162) + cos(t163);
	t158 = t174 * t165 + t170 * t179;
	t157 = -t177 * t165 + t170 * t180;
	t156 = t160 * t176 + t161 * t173 + pkin(1);
	t155 = t168 * t171 - t178 * t170;
	t1 = [-t157 * t175 + t174 * t183, t157 * t172 + t174 * t182, t179 + t174 * t159 / 0.2e1, t155 * t174 + t156 * t177 + 0; t158 * t175 - t172 * t181, -t158 * t172 - t175 * t181, t180 - t177 * t159 / 0.2e1, -t155 * t177 + t156 * t174 + 0; t164 * t182 + t170 * t172, -t164 * t183 + t170 * t175, -sin(t163) / 0.2e1 - sin(t162) / 0.2e1, t178 * t168 + t170 * t171 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:14:34
	% EndTime: 2020-11-04 22:14:34
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (102->44), mult. (170->73), div. (0->0), fcn. (204->14), ass. (0->39)
	t191 = qJ(2) + pkin(11);
	t189 = sin(t191);
	t195 = cos(pkin(6));
	t221 = t189 * t195;
	t193 = sin(pkin(6));
	t197 = sin(qJ(4));
	t220 = t193 * t197;
	t200 = cos(qJ(5));
	t219 = t193 * t200;
	t196 = sin(qJ(5));
	t199 = sin(qJ(1));
	t218 = t196 * t199;
	t203 = cos(qJ(1));
	t217 = t196 * t203;
	t201 = cos(qJ(4));
	t216 = t199 * t201;
	t215 = t200 * t199;
	t214 = t200 * t203;
	t213 = t201 * t203;
	t212 = t196 * t220;
	t211 = t197 * t219;
	t210 = t201 * t215;
	t209 = t196 * t216;
	t208 = t196 * t213;
	t207 = t200 * t213;
	t192 = sin(pkin(11));
	t194 = cos(pkin(11));
	t205 = pkin(4) * t201 + pkin(10) * t197 + pkin(3);
	t186 = t192 * pkin(9) + t205 * t194 + pkin(2);
	t187 = t194 * pkin(9) - t205 * t192;
	t198 = sin(qJ(2));
	t202 = cos(qJ(2));
	t206 = t186 * t198 - t187 * t202;
	t204 = t197 * pkin(4) - pkin(10) * t201 + pkin(8) + qJ(3);
	t190 = cos(t191);
	t188 = t193 * t189 * t201 + t195 * t197;
	t185 = t186 * t202 + t187 * t198 + pkin(1);
	t184 = t193 * t204 - t206 * t195;
	t1 = [(t195 * t218 + t207) * t190 + (-t195 * t210 + t217) * t189 + t199 * t211, (t195 * t215 - t208) * t190 + (t195 * t209 + t214) * t189 - t199 * t212, -(-t203 * t190 + t199 * t221) * t197 - t193 * t216, t184 * t199 + t185 * t203 + 0; (-t195 * t217 + t210) * t190 + (t195 * t207 + t218) * t189 - t203 * t211, (-t195 * t214 - t209) * t190 + (-t195 * t208 + t215) * t189 + t203 * t212, (t199 * t190 + t203 * t221) * t197 + t193 * t213, -t184 * t203 + t185 * t199 + 0; -t193 * t190 * t196 + t188 * t200, -t188 * t196 - t190 * t219, t189 * t220 - t195 * t201, t206 * t193 + t204 * t195 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:14:34
	% EndTime: 2020-11-04 22:14:34
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (148->48), mult. (242->77), div. (0->0), fcn. (276->14), ass. (0->41)
	t236 = sin(qJ(5));
	t240 = cos(qJ(5));
	t227 = pkin(5) * t240 + qJ(6) * t236 + pkin(4);
	t237 = sin(qJ(4));
	t241 = cos(qJ(4));
	t265 = -pkin(10) * t241 + t227 * t237 + pkin(8) + qJ(3);
	t231 = qJ(2) + pkin(11);
	t228 = sin(t231);
	t235 = cos(pkin(6));
	t263 = t228 * t235;
	t233 = sin(pkin(6));
	t262 = t233 * t237;
	t261 = t233 * t240;
	t239 = sin(qJ(1));
	t260 = t236 * t239;
	t243 = cos(qJ(1));
	t259 = t236 * t243;
	t258 = t239 * t241;
	t257 = t240 * t239;
	t256 = t240 * t243;
	t255 = t241 * t243;
	t254 = t236 * t262;
	t253 = t237 * t261;
	t252 = t241 * t257;
	t251 = t236 * t258;
	t250 = t236 * t255;
	t249 = t240 * t255;
	t232 = sin(pkin(11));
	t234 = cos(pkin(11));
	t244 = pkin(5) * t236 - qJ(6) * t240 + pkin(9);
	t245 = pkin(10) * t237 + t227 * t241 + pkin(3);
	t224 = t244 * t232 + t245 * t234 + pkin(2);
	t225 = t245 * t232 - t244 * t234;
	t238 = sin(qJ(2));
	t242 = cos(qJ(2));
	t246 = t224 * t238 + t225 * t242;
	t229 = cos(t231);
	t226 = t233 * t228 * t241 + t235 * t237;
	t223 = t224 * t242 - t225 * t238 + pkin(1);
	t222 = -t233 * t265 + t246 * t235;
	t1 = [(t235 * t260 + t249) * t229 + (-t235 * t252 + t259) * t228 + t239 * t253, -(-t243 * t229 + t239 * t263) * t237 - t233 * t258, (-t235 * t257 + t250) * t229 + (-t235 * t251 - t256) * t228 + t239 * t254, -t222 * t239 + t223 * t243 + 0; (-t235 * t259 + t252) * t229 + (t235 * t249 + t260) * t228 - t243 * t253, (t239 * t229 + t243 * t263) * t237 + t233 * t255, (t235 * t256 + t251) * t229 + (t235 * t250 - t257) * t228 - t243 * t254, t222 * t243 + t223 * t239 + 0; -t233 * t229 * t236 + t226 * t240, t228 * t262 - t235 * t241, t226 * t236 + t229 * t261, t246 * t233 + t265 * t235 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end