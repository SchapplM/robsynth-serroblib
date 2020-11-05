% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRPRPR3 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:59
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6PRPRPR3_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR3_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRPR3_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR3_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:59:34
	% EndTime: 2020-11-04 20:59:34
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:59:34
	% EndTime: 2020-11-04 20:59:34
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t121 = cos(pkin(10));
	t120 = sin(pkin(10));
	t1 = [t121, -t120, 0, 0; t120, t121, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:59:34
	% EndTime: 2020-11-04 20:59:34
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->19), div. (0->0), fcn. (36->6), ass. (0->11)
	t122 = sin(pkin(10));
	t123 = sin(pkin(6));
	t131 = t122 * t123;
	t124 = cos(pkin(10));
	t130 = t124 * t123;
	t125 = cos(pkin(6));
	t126 = sin(qJ(2));
	t129 = t125 * t126;
	t127 = cos(qJ(2));
	t128 = t125 * t127;
	t1 = [-t122 * t129 + t124 * t127, -t122 * t128 - t124 * t126, t131, t124 * pkin(1) + pkin(7) * t131 + 0; t122 * t127 + t124 * t129, -t122 * t126 + t124 * t128, -t130, t122 * pkin(1) - pkin(7) * t130 + 0; t123 * t126, t123 * t127, t125, t125 * pkin(7) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:59:34
	% EndTime: 2020-11-04 20:59:34
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (53->22), mult. (37->26), div. (0->0), fcn. (44->12), ass. (0->20)
	t151 = pkin(2) * sin(qJ(2));
	t144 = qJ(2) + pkin(11);
	t149 = pkin(7) + qJ(3);
	t148 = cos(pkin(6));
	t147 = cos(pkin(10));
	t146 = sin(pkin(6));
	t145 = sin(pkin(10));
	t143 = cos(t144);
	t142 = sin(t144);
	t141 = pkin(6) - t144;
	t140 = pkin(6) + t144;
	t139 = cos(qJ(2)) * pkin(2) + pkin(1);
	t138 = cos(t140);
	t137 = sin(t141);
	t136 = cos(t141) / 0.2e1;
	t135 = sin(t140) / 0.2e1;
	t134 = -t146 * t149 + t148 * t151;
	t133 = t138 / 0.2e1 + t136;
	t132 = t135 - t137 / 0.2e1;
	t1 = [-t145 * t132 + t147 * t143, -t145 * t133 - t147 * t142, t145 * t146, -t145 * t134 + t147 * t139 + 0; t147 * t132 + t145 * t143, t147 * t133 - t145 * t142, -t147 * t146, t147 * t134 + t145 * t139 + 0; t136 - t138 / 0.2e1, t137 / 0.2e1 + t135, t148, t146 * t151 + t148 * t149 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:59:34
	% EndTime: 2020-11-04 20:59:34
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (73->35), mult. (85->53), div. (0->0), fcn. (105->16), ass. (0->27)
	t163 = sin(pkin(10));
	t167 = cos(pkin(6));
	t177 = t163 * t167;
	t164 = sin(pkin(6));
	t168 = pkin(7) + qJ(3);
	t176 = t164 * t168;
	t169 = sin(qJ(4));
	t175 = t164 * t169;
	t171 = cos(qJ(4));
	t174 = t164 * t171;
	t166 = cos(pkin(10));
	t173 = t166 * t167;
	t161 = qJ(2) + pkin(11);
	t172 = cos(qJ(2));
	t170 = sin(qJ(2));
	t165 = cos(pkin(11));
	t162 = sin(pkin(11));
	t160 = cos(t161);
	t159 = sin(t161);
	t158 = pkin(6) - t161;
	t157 = pkin(6) + t161;
	t156 = -t162 * pkin(3) + t165 * pkin(8);
	t155 = t165 * pkin(3) + t162 * pkin(8) + pkin(2);
	t154 = cos(t157) + cos(t158);
	t153 = t159 * t173 + t163 * t160;
	t152 = t159 * t177 - t166 * t160;
	t1 = [-t152 * t171 + t163 * t175, t152 * t169 + t163 * t174, t166 * t159 + t163 * t154 / 0.2e1, (t166 * t155 + t156 * t177) * t172 + (-t155 * t177 + t166 * t156) * t170 + t163 * t176 + t166 * pkin(1) + 0; t153 * t171 - t166 * t175, -t153 * t169 - t166 * t174, t163 * t159 - t166 * t154 / 0.2e1, (t163 * t155 - t156 * t173) * t172 + (t155 * t173 + t163 * t156) * t170 - t166 * t176 + t163 * pkin(1) + 0; t159 * t174 + t167 * t169, -t159 * t175 + t167 * t171, -sin(t158) / 0.2e1 - sin(t157) / 0.2e1, t167 * t168 + qJ(1) + 0 + (t155 * t170 - t156 * t172) * t164; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:59:34
	% EndTime: 2020-11-04 20:59:34
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (99->45), mult. (121->59), div. (0->0), fcn. (151->16), ass. (0->33)
	t195 = sin(pkin(10));
	t199 = cos(pkin(6));
	t209 = t195 * t199;
	t196 = sin(pkin(6));
	t200 = pkin(7) + qJ(3);
	t208 = t196 * t200;
	t201 = sin(qJ(4));
	t207 = t196 * t201;
	t203 = cos(qJ(4));
	t206 = t196 * t203;
	t198 = cos(pkin(10));
	t205 = t198 * t199;
	t193 = qJ(2) + pkin(11);
	t204 = cos(qJ(2));
	t202 = sin(qJ(2));
	t197 = cos(pkin(11));
	t194 = sin(pkin(11));
	t192 = cos(t193);
	t191 = sin(t193);
	t190 = pkin(6) - t193;
	t189 = pkin(6) + t193;
	t188 = -t194 * pkin(3) + t197 * pkin(8);
	t187 = t197 * pkin(3) + t194 * pkin(8) + pkin(2);
	t186 = cos(t189) + cos(t190);
	t185 = t191 * t206 + t199 * t201;
	t184 = t191 * t207 - t199 * t203;
	t183 = t191 * t205 + t195 * t192;
	t182 = t191 * t209 - t198 * t192;
	t181 = -t182 * t203 + t195 * t207;
	t180 = t183 * t203 - t198 * t207;
	t179 = t183 * t201 + t198 * t206;
	t178 = t182 * t201 + t195 * t206;
	t1 = [t198 * t191 + t195 * t186 / 0.2e1, -t181, -t178, t181 * pkin(4) - t178 * qJ(5) + (t198 * t187 + t188 * t209) * t204 + (-t187 * t209 + t198 * t188) * t202 + t195 * t208 + t198 * pkin(1) + 0; t195 * t191 - t198 * t186 / 0.2e1, -t180, t179, t180 * pkin(4) + t179 * qJ(5) + (t195 * t187 - t188 * t205) * t204 + (t187 * t205 + t195 * t188) * t202 - t198 * t208 + t195 * pkin(1) + 0; -sin(t190) / 0.2e1 - sin(t189) / 0.2e1, -t185, t184, t185 * pkin(4) + t184 * qJ(5) + t199 * t200 + qJ(1) + 0 + (t187 * t202 - t188 * t204) * t196; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:59:34
	% EndTime: 2020-11-04 20:59:35
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (123->65), mult. (194->99), div. (0->0), fcn. (231->16), ass. (0->44)
	t223 = sin(pkin(6));
	t229 = sin(qJ(4));
	t232 = cos(qJ(4));
	t234 = pkin(9) + pkin(4);
	t236 = qJ(5) * t232 - t229 * t234 - pkin(7) - qJ(3);
	t260 = t236 * t223;
	t237 = qJ(5) * t229 + t234 * t232;
	t221 = sin(pkin(11));
	t259 = t221 * pkin(3);
	t224 = cos(pkin(11));
	t258 = t224 * pkin(3) + pkin(2);
	t220 = qJ(2) + pkin(11);
	t218 = cos(t220);
	t256 = t218 * t223;
	t222 = sin(pkin(10));
	t255 = t222 * t221;
	t228 = sin(qJ(6));
	t254 = t222 * t228;
	t231 = cos(qJ(6));
	t253 = t222 * t231;
	t252 = t223 * t229;
	t251 = t224 * t222;
	t225 = cos(pkin(10));
	t226 = cos(pkin(6));
	t250 = t225 * t226;
	t249 = t226 * t232;
	t248 = t228 * t225;
	t247 = t228 * t229;
	t246 = t231 * t225;
	t245 = t232 * t223;
	t243 = t222 * t247;
	t242 = t229 * t253;
	t241 = t225 * t247;
	t240 = t226 * t246;
	t239 = t228 * t245;
	t238 = t231 * t245;
	t235 = pkin(5) + pkin(8);
	t233 = cos(qJ(2));
	t230 = sin(qJ(2));
	t217 = sin(t220);
	t216 = t235 * t224 - t259;
	t215 = t235 * t221 + t258;
	t210 = t217 * t252 - t249;
	t1 = [(t226 * t253 + t241) * t218 + (-t226 * t243 + t246) * t217 - t222 * t239, (-t226 * t254 + t229 * t246) * t218 + (-t226 * t242 - t248) * t217 - t222 * t238, t232 * t225 * t218 + (-t217 * t249 + t252) * t222, (t215 * t225 - t237 * (-t225 * t224 + t226 * t255)) * t233 + (t225 * t216 - t237 * (t225 * t221 + t226 * t251)) * t230 + t225 * pkin(1) + 0 + ((-t215 * t230 + t216 * t233) * t226 - t260) * t222; (-t240 + t243) * t218 + (t226 * t241 + t253) * t217 + t225 * t239, (t226 * t248 + t242) * t218 + (t229 * t240 - t254) * t217 + t225 * t238, (t217 * t250 + t222 * t218) * t232 - t225 * t252, (-t216 * t250 + t215 * t222 + t237 * (t221 * t250 + t251)) * t233 + (t215 * t250 + t222 * t216 + t237 * (t224 * t250 - t255)) * t230 + t222 * pkin(1) + 0 + t225 * t260; t210 * t228 - t231 * t256, t210 * t231 + t228 * t256, t217 * t245 + t226 * t229, qJ(1) + 0 + (-sin(pkin(6) - t220) / 0.2e1 - sin(pkin(6) + t220) / 0.2e1) * pkin(5) - t236 * t226 + (t237 * t217 + (t221 * pkin(8) + t258) * t230 - (t224 * pkin(8) - t259) * t233) * t223; 0, 0, 0, 1;];
	Tc_mdh = t1;
end