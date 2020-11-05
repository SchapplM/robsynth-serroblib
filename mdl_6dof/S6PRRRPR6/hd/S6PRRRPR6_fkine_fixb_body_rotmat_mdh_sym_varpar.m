% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRRRPR6 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:16
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6PRRRPR6_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR6_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRPR6_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR6_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:16:37
	% EndTime: 2020-11-04 21:16:37
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:16:37
	% EndTime: 2020-11-04 21:16:37
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t137 = cos(pkin(11));
	t136 = sin(pkin(11));
	t1 = [t137, -t136, 0, 0; t136, t137, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:16:37
	% EndTime: 2020-11-04 21:16:37
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->19), div. (0->0), fcn. (36->6), ass. (0->11)
	t138 = sin(pkin(11));
	t139 = sin(pkin(6));
	t147 = t138 * t139;
	t140 = cos(pkin(11));
	t146 = t140 * t139;
	t141 = cos(pkin(6));
	t142 = sin(qJ(2));
	t145 = t141 * t142;
	t143 = cos(qJ(2));
	t144 = t141 * t143;
	t1 = [-t138 * t145 + t140 * t143, -t138 * t144 - t140 * t142, t147, t140 * pkin(1) + pkin(7) * t147 + 0; t138 * t143 + t140 * t145, -t138 * t142 + t140 * t144, -t146, t138 * pkin(1) - pkin(7) * t146 + 0; t139 * t142, t139 * t143, t141, t141 * pkin(7) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:16:37
	% EndTime: 2020-11-04 21:16:37
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (29->27), mult. (64->48), div. (0->0), fcn. (85->8), ass. (0->18)
	t151 = sin(pkin(6));
	t164 = t151 * pkin(7);
	t150 = sin(pkin(11));
	t153 = cos(pkin(6));
	t163 = t150 * t153;
	t154 = sin(qJ(3));
	t162 = t151 * t154;
	t156 = cos(qJ(3));
	t161 = t151 * t156;
	t152 = cos(pkin(11));
	t160 = t152 * t153;
	t155 = sin(qJ(2));
	t159 = t153 * t155;
	t157 = cos(qJ(2));
	t158 = t153 * t157;
	t149 = -t150 * t159 + t152 * t157;
	t148 = t150 * t157 + t152 * t159;
	t1 = [t149 * t156 + t150 * t162, -t149 * t154 + t150 * t161, t150 * t158 + t152 * t155, (t152 * pkin(2) + pkin(8) * t163) * t157 + (-pkin(2) * t163 + t152 * pkin(8)) * t155 + t150 * t164 + t152 * pkin(1) + 0; t148 * t156 - t152 * t162, -t148 * t154 - t152 * t161, t150 * t155 - t152 * t158, (t150 * pkin(2) - pkin(8) * t160) * t157 + (pkin(2) * t160 + t150 * pkin(8)) * t155 - t152 * t164 + t150 * pkin(1) + 0; t153 * t154 + t155 * t161, t153 * t156 - t155 * t162, -t151 * t157, t153 * pkin(7) + qJ(1) + 0 + (pkin(2) * t155 - pkin(8) * t157) * t151; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:16:37
	% EndTime: 2020-11-04 21:16:37
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (57->32), mult. (136->56), div. (0->0), fcn. (170->10), ass. (0->29)
	t172 = sin(pkin(6));
	t174 = cos(pkin(6));
	t177 = sin(qJ(2));
	t180 = cos(qJ(2));
	t176 = sin(qJ(3));
	t179 = cos(qJ(3));
	t184 = t179 * pkin(3) + t176 * pkin(9) + pkin(2);
	t182 = -pkin(8) * t180 + t184 * t177;
	t183 = t176 * pkin(3) - t179 * pkin(9) + pkin(7);
	t195 = t183 * t172 - t182 * t174;
	t191 = t172 * t179;
	t190 = t172 * t180;
	t189 = t174 * t177;
	t188 = t174 * t180;
	t187 = t176 * t172;
	t186 = t177 * t179;
	t185 = t179 * t180;
	t181 = pkin(8) * t177 + t184 * t180 + pkin(1);
	t178 = cos(qJ(4));
	t175 = sin(qJ(4));
	t173 = cos(pkin(11));
	t171 = sin(pkin(11));
	t170 = t174 * t186 - t187;
	t169 = t172 * t186 + t174 * t176;
	t168 = t171 * t188 + t173 * t177;
	t167 = t171 * t177 - t173 * t188;
	t166 = -t171 * t170 + t173 * t185;
	t165 = t173 * t170 + t171 * t185;
	t1 = [t166 * t178 + t168 * t175, -t166 * t175 + t168 * t178, -(t171 * t189 - t173 * t180) * t176 - t171 * t191, t195 * t171 + t181 * t173 + 0; t165 * t178 + t167 * t175, -t165 * t175 + t167 * t178, (t171 * t180 + t173 * t189) * t176 + t173 * t191, t181 * t171 - t195 * t173 + 0; t169 * t178 - t175 * t190, -t169 * t175 - t178 * t190, -t174 * t179 + t177 * t187, t182 * t172 + t183 * t174 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:16:37
	% EndTime: 2020-11-04 21:16:37
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (83->36), mult. (180->60), div. (0->0), fcn. (214->10), ass. (0->31)
	t206 = sin(pkin(6));
	t208 = cos(pkin(6));
	t211 = sin(qJ(2));
	t214 = cos(qJ(2));
	t209 = sin(qJ(4));
	t212 = cos(qJ(4));
	t217 = t209 * pkin(4) - qJ(5) * t212 + pkin(8);
	t210 = sin(qJ(3));
	t213 = cos(qJ(3));
	t230 = pkin(4) * t212 + qJ(5) * t209 + pkin(3);
	t218 = t210 * pkin(9) + t213 * t230 + pkin(2);
	t216 = t218 * t211 - t217 * t214;
	t229 = -t213 * pkin(9) + t210 * t230 + pkin(7);
	t232 = t229 * t206 - t216 * t208;
	t227 = t206 * t213;
	t226 = t206 * t214;
	t225 = t208 * t211;
	t224 = t208 * t214;
	t223 = t210 * t206;
	t222 = t211 * t213;
	t221 = t213 * t214;
	t215 = t217 * t211 + t218 * t214 + pkin(1);
	t207 = cos(pkin(11));
	t205 = sin(pkin(11));
	t201 = t208 * t222 - t223;
	t200 = t206 * t222 + t208 * t210;
	t199 = t205 * t224 + t207 * t211;
	t198 = t205 * t211 - t207 * t224;
	t197 = -t205 * t201 + t207 * t221;
	t196 = t207 * t201 + t205 * t221;
	t1 = [t197 * t212 + t199 * t209, -(t205 * t225 - t207 * t214) * t210 - t205 * t227, t197 * t209 - t199 * t212, t232 * t205 + t215 * t207 + 0; t196 * t212 + t198 * t209, (t205 * t214 + t207 * t225) * t210 + t207 * t227, t196 * t209 - t198 * t212, t215 * t205 - t232 * t207 + 0; t200 * t212 - t209 * t226, -t208 * t213 + t211 * t223, t200 * t209 + t212 * t226, t216 * t206 + t229 * t208 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:16:37
	% EndTime: 2020-11-04 21:16:38
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (130->44), mult. (250->72), div. (0->0), fcn. (310->12), ass. (0->41)
	t247 = sin(pkin(6));
	t249 = cos(pkin(6));
	t253 = sin(qJ(2));
	t257 = cos(qJ(2));
	t251 = sin(qJ(4));
	t255 = cos(qJ(4));
	t259 = pkin(4) + pkin(5);
	t263 = qJ(5) * t255 - t259 * t251 - pkin(8);
	t252 = sin(qJ(3));
	t256 = cos(qJ(3));
	t258 = pkin(10) - pkin(9);
	t278 = qJ(5) * t251 + t259 * t255 + pkin(3);
	t264 = -t258 * t252 + t256 * t278 + pkin(2);
	t261 = t264 * t253 + t263 * t257;
	t277 = t252 * t278 + t258 * t256 + pkin(7);
	t280 = t277 * t247 - t261 * t249;
	t276 = t247 * t256;
	t275 = t247 * t257;
	t274 = t249 * t253;
	t273 = t249 * t257;
	t272 = t252 * t247;
	t271 = t253 * t256;
	t270 = t256 * t257;
	t246 = sin(pkin(11));
	t248 = cos(pkin(11));
	t262 = -t249 * t271 + t272;
	t235 = t246 * t270 - t248 * t262;
	t240 = t246 * t253 - t248 * t273;
	t250 = sin(qJ(6));
	t254 = cos(qJ(6));
	t266 = t235 * t254 - t250 * t240;
	t236 = t246 * t262 + t248 * t270;
	t241 = t246 * t273 + t248 * t253;
	t265 = t250 * t236 + t241 * t254;
	t260 = -t263 * t253 + t264 * t257 + pkin(1);
	t242 = t247 * t271 + t249 * t252;
	t238 = t242 * t254 + t250 * t275;
	t237 = -t250 * t242 + t254 * t275;
	t234 = -t250 * t235 - t240 * t254;
	t233 = t236 * t254 - t250 * t241;
	t1 = [t233 * t255 + t251 * t265, t251 * t233 - t265 * t255, (t246 * t274 - t248 * t257) * t252 + t246 * t276, t280 * t246 + t260 * t248 + 0; -t251 * t234 + t266 * t255, t234 * t255 + t251 * t266, -(t246 * t257 + t248 * t274) * t252 - t248 * t276, t260 * t246 - t280 * t248 + 0; -t251 * t237 + t238 * t255, t237 * t255 + t238 * t251, t249 * t256 - t253 * t272, t261 * t247 + t277 * t249 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end