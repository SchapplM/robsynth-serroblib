% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRRRR12 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d1,d2,d3,d4,d5]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:52
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RRRRR12_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR12_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRR12_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RRRRR12_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:52:08
	% EndTime: 2020-11-04 20:52:08
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:52:08
	% EndTime: 2020-11-04 20:52:08
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t137 = cos(qJ(1));
	t136 = sin(qJ(1));
	t1 = [t137, -t136, 0, 0; t136, t137, 0, 0; 0, 0, 1, pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:52:08
	% EndTime: 2020-11-04 20:52:08
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->17), div. (0->0), fcn. (36->6), ass. (0->13)
	t138 = sin(pkin(5));
	t141 = sin(qJ(1));
	t149 = t141 * t138;
	t140 = sin(qJ(2));
	t148 = t141 * t140;
	t142 = cos(qJ(2));
	t147 = t141 * t142;
	t143 = cos(qJ(1));
	t146 = t143 * t138;
	t145 = t143 * t140;
	t144 = t143 * t142;
	t139 = cos(pkin(5));
	t1 = [-t139 * t148 + t144, -t139 * t147 - t145, t149, t143 * pkin(1) + pkin(8) * t149 + 0; t139 * t145 + t147, t139 * t144 - t148, -t146, t141 * pkin(1) - pkin(8) * t146 + 0; t138 * t140, t138 * t142, t139, t139 * pkin(8) + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:52:08
	% EndTime: 2020-11-04 20:52:08
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (44->32), mult. (104->58), div. (0->0), fcn. (138->10), ass. (0->29)
	t156 = sin(pkin(6));
	t159 = cos(pkin(5));
	t178 = t156 * t159;
	t164 = cos(qJ(2));
	t177 = t156 * t164;
	t157 = sin(pkin(5));
	t158 = cos(pkin(6));
	t176 = t157 * t158;
	t175 = t159 * t164;
	t160 = sin(qJ(3));
	t161 = sin(qJ(2));
	t174 = t160 * t161;
	t173 = t160 * t164;
	t163 = cos(qJ(3));
	t172 = t161 * t163;
	t162 = sin(qJ(1));
	t171 = t162 * t161;
	t170 = t164 * t163;
	t165 = cos(qJ(1));
	t169 = t165 * t161;
	t168 = t165 * t164;
	t167 = -pkin(2) * t161 + pkin(9) * t177;
	t151 = -t156 * t157 + t158 * t175;
	t166 = t151 * t160 + t159 * t172;
	t155 = t158 * pkin(9) + pkin(8);
	t153 = t156 * t161 * pkin(9) + pkin(2) * t164 + pkin(1);
	t152 = -t158 * t174 + t170;
	t150 = t157 * t155 + t167 * t159;
	t1 = [t165 * t152 - t166 * t162, (-t151 * t162 - t158 * t169) * t163 - (-t159 * t171 + t168) * t160, (t162 * t175 + t169) * t156 + t162 * t176, t150 * t162 + t153 * t165 + 0; t162 * t152 + t166 * t165, (t151 * t163 - t159 * t174) * t165 - t162 * (t158 * t172 + t173), -(t159 * t168 - t171) * t156 - t165 * t176, -t150 * t165 + t153 * t162 + 0; t160 * t178 + (t158 * t173 + t172) * t157, t163 * t178 + (t158 * t170 - t174) * t157, -t157 * t177 + t159 * t158, t155 * t159 - t167 * t157 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:52:08
	% EndTime: 2020-11-04 20:52:08
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (91->44), mult. (209->80), div. (0->0), fcn. (264->12), ass. (0->38)
	t194 = sin(pkin(6));
	t196 = cos(pkin(6));
	t199 = sin(qJ(3));
	t203 = cos(qJ(3));
	t209 = pkin(3) * t199 - pkin(10) * t203;
	t184 = t194 * pkin(9) - t209 * t196;
	t189 = t203 * pkin(3) + pkin(10) * t199 + pkin(2);
	t200 = sin(qJ(2));
	t204 = cos(qJ(2));
	t221 = -t184 * t204 + t189 * t200;
	t195 = sin(pkin(5));
	t218 = t194 * t195;
	t217 = t194 * t199;
	t216 = t194 * t200;
	t197 = cos(pkin(5));
	t215 = t197 * t204;
	t214 = t199 * t200;
	t213 = t199 * t204;
	t212 = t200 * t203;
	t205 = cos(qJ(1));
	t211 = t200 * t205;
	t210 = t204 * t203;
	t207 = t196 * t213 + t212;
	t181 = -t195 * t217 + t207 * t197;
	t185 = t194 * t215 + t195 * t196;
	t198 = sin(qJ(4));
	t202 = cos(qJ(4));
	t208 = t181 * t198 + t202 * t185;
	t206 = t196 * pkin(9) + t209 * t194 + pkin(8);
	t201 = sin(qJ(1));
	t188 = t196 * t214 - t210;
	t187 = t196 * t215 - t218;
	t186 = -t197 * t196 + t204 * t218;
	t183 = t188 * t198 + t202 * t216;
	t182 = t207 * t195 + t197 * t217;
	t180 = t184 * t200 + t189 * t204 + pkin(1);
	t179 = t195 * t206 - t221 * t197;
	t1 = [(-t181 * t201 - t205 * t188) * t202 + (t185 * t201 + t194 * t211) * t198, t183 * t205 + t208 * t201, (t187 * t201 + t196 * t211) * t203 + (-t201 * t197 * t200 + t205 * t204) * t199, t179 * t201 + t180 * t205 + 0; (t181 * t202 - t198 * t185) * t205 + (-t188 * t202 + t198 * t216) * t201, t183 * t201 - t208 * t205, (-t187 * t203 + t197 * t214) * t205 + t201 * (t196 * t212 + t213), -t179 * t205 + t180 * t201 + 0; t182 * t202 - t198 * t186, -t182 * t198 - t202 * t186, -t197 * t194 * t203 + (-t196 * t210 + t214) * t195, t221 * t195 + t206 * t197 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:52:08
	% EndTime: 2020-11-04 20:52:09
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (167->55), mult. (361->96), div. (0->0), fcn. (450->14), ass. (0->48)
	t242 = sin(pkin(6));
	t248 = sin(qJ(3));
	t253 = cos(qJ(3));
	t247 = sin(qJ(4));
	t252 = cos(qJ(4));
	t278 = t252 * pkin(4) + t247 * pkin(11) + pkin(3);
	t279 = -pkin(10) * t253 + t248 * t278;
	t281 = -t279 * t242 - pkin(8);
	t244 = cos(pkin(6));
	t260 = t247 * pkin(4) - t252 * pkin(11) + pkin(9);
	t227 = t242 * t260 - t244 * t279;
	t235 = pkin(10) * t248 + t253 * t278 + pkin(2);
	t249 = sin(qJ(2));
	t254 = cos(qJ(2));
	t280 = t227 * t254 - t235 * t249;
	t274 = t242 * t248;
	t273 = t242 * t253;
	t272 = t242 * t254;
	t243 = sin(pkin(5));
	t271 = t243 * t244;
	t270 = t248 * t249;
	t269 = t248 * t252;
	t268 = t248 * t254;
	t267 = t249 * t253;
	t266 = t254 * t253;
	t245 = cos(pkin(5));
	t258 = t242 * t269 + t244 * t247;
	t237 = -t242 * t247 + t244 * t269;
	t259 = t237 * t254 + t252 * t267;
	t225 = -t243 * t258 + t259 * t245;
	t257 = t244 * t266 - t270;
	t229 = -t243 * t273 + t257 * t245;
	t246 = sin(qJ(5));
	t251 = cos(qJ(5));
	t262 = t225 * t246 + t251 * t229;
	t256 = t244 * t268 + t267;
	t261 = (-t243 * t274 + t256 * t245) * t247 + t252 * (t245 * t272 + t271);
	t255 = cos(qJ(1));
	t250 = sin(qJ(1));
	t238 = t244 * t267 + t268;
	t232 = (t244 * t270 - t266) * t247 + t242 * t252 * t249;
	t231 = -t249 * t237 + t252 * t266;
	t230 = t257 * t243 + t245 * t273;
	t226 = -t231 * t246 + t251 * t238;
	t224 = -t259 * t243 - t245 * t258;
	t223 = t227 * t249 + t235 * t254 + pkin(1);
	t222 = -t281 * t243 + t280 * t245 + t260 * t271;
	t1 = [(-t225 * t250 + t255 * t231) * t251 + (t229 * t250 + t255 * t238) * t246, t226 * t255 + t262 * t250, -t232 * t255 - t261 * t250, t222 * t250 + t223 * t255 + 0; (t225 * t251 - t246 * t229) * t255 + t250 * (t231 * t251 + t246 * t238), t226 * t250 - t262 * t255, -t232 * t250 + t261 * t255, -t222 * t255 + t223 * t250 + 0; -t224 * t251 - t246 * t230, t224 * t246 - t251 * t230, (t256 * t243 + t245 * t274) * t247 + t252 * (t243 * t272 - t245 * t244), pkin(7) + 0 - t280 * t243 + (t260 * t244 - t281) * t245; 0, 0, 0, 1;];
	Tc_mdh = t1;
end