% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRPRPR1 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:58
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6PRPRPR1_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR1_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRPR1_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR1_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [12x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:58:39
	% EndTime: 2020-11-04 20:58:39
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:58:39
	% EndTime: 2020-11-04 20:58:39
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t131 = cos(pkin(10));
	t130 = sin(pkin(10));
	t1 = [t131, -t130, 0, 0; t130, t131, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:58:39
	% EndTime: 2020-11-04 20:58:39
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->19), div. (0->0), fcn. (36->6), ass. (0->11)
	t132 = sin(pkin(10));
	t133 = sin(pkin(6));
	t141 = t132 * t133;
	t134 = cos(pkin(10));
	t140 = t134 * t133;
	t135 = cos(pkin(6));
	t136 = sin(qJ(2));
	t139 = t135 * t136;
	t137 = cos(qJ(2));
	t138 = t135 * t137;
	t1 = [-t132 * t139 + t134 * t137, -t132 * t138 - t134 * t136, t141, t134 * pkin(1) + pkin(7) * t141 + 0; t132 * t137 + t134 * t139, -t132 * t136 + t134 * t138, -t140, t132 * pkin(1) - pkin(7) * t140 + 0; t133 * t136, t133 * t137, t135, t135 * pkin(7) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:58:39
	% EndTime: 2020-11-04 20:58:39
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (53->22), mult. (37->26), div. (0->0), fcn. (44->12), ass. (0->20)
	t161 = pkin(2) * sin(qJ(2));
	t154 = qJ(2) + pkin(11);
	t159 = pkin(7) + qJ(3);
	t158 = cos(pkin(6));
	t157 = cos(pkin(10));
	t156 = sin(pkin(6));
	t155 = sin(pkin(10));
	t153 = cos(t154);
	t152 = sin(t154);
	t151 = pkin(6) - t154;
	t150 = pkin(6) + t154;
	t149 = cos(qJ(2)) * pkin(2) + pkin(1);
	t148 = cos(t150);
	t147 = sin(t151);
	t146 = cos(t151) / 0.2e1;
	t145 = sin(t150) / 0.2e1;
	t144 = -t156 * t159 + t158 * t161;
	t143 = t148 / 0.2e1 + t146;
	t142 = t145 - t147 / 0.2e1;
	t1 = [-t155 * t142 + t157 * t153, -t155 * t143 - t157 * t152, t155 * t156, -t155 * t144 + t157 * t149 + 0; t157 * t142 + t155 * t153, t157 * t143 - t155 * t152, -t157 * t156, t157 * t144 + t155 * t149 + 0; t146 - t148 / 0.2e1, t147 / 0.2e1 + t145, t158, t156 * t161 + t158 * t159 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:58:39
	% EndTime: 2020-11-04 20:58:39
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (73->35), mult. (85->53), div. (0->0), fcn. (105->16), ass. (0->27)
	t173 = sin(pkin(10));
	t177 = cos(pkin(6));
	t187 = t173 * t177;
	t174 = sin(pkin(6));
	t178 = pkin(7) + qJ(3);
	t186 = t174 * t178;
	t179 = sin(qJ(4));
	t185 = t174 * t179;
	t181 = cos(qJ(4));
	t184 = t174 * t181;
	t176 = cos(pkin(10));
	t183 = t176 * t177;
	t171 = qJ(2) + pkin(11);
	t182 = cos(qJ(2));
	t180 = sin(qJ(2));
	t175 = cos(pkin(11));
	t172 = sin(pkin(11));
	t170 = cos(t171);
	t169 = sin(t171);
	t168 = pkin(6) - t171;
	t167 = pkin(6) + t171;
	t166 = -t172 * pkin(3) + t175 * pkin(8);
	t165 = t175 * pkin(3) + t172 * pkin(8) + pkin(2);
	t164 = cos(t167) + cos(t168);
	t163 = t169 * t183 + t173 * t170;
	t162 = t169 * t187 - t176 * t170;
	t1 = [-t162 * t181 + t173 * t185, t162 * t179 + t173 * t184, t176 * t169 + t173 * t164 / 0.2e1, (t176 * t165 + t166 * t187) * t182 + (-t165 * t187 + t176 * t166) * t180 + t173 * t186 + t176 * pkin(1) + 0; t163 * t181 - t176 * t185, -t163 * t179 - t176 * t184, t173 * t169 - t176 * t164 / 0.2e1, (t173 * t165 - t166 * t183) * t182 + (t165 * t183 + t173 * t166) * t180 - t176 * t186 + t173 * pkin(1) + 0; t169 * t184 + t177 * t179, -t169 * t185 + t177 * t181, -sin(t168) / 0.2e1 - sin(t167) / 0.2e1, t177 * t178 + qJ(1) + 0 + (t165 * t180 - t166 * t182) * t174; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:58:39
	% EndTime: 2020-11-04 20:58:39
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (142->38), mult. (98->46), div. (0->0), fcn. (108->16), ass. (0->33)
	t226 = pkin(7) + qJ(3) + pkin(4) * sin(qJ(4));
	t210 = qJ(2) + pkin(11);
	t203 = pkin(6) + t210;
	t225 = sin(t203) / 0.2e1;
	t204 = pkin(6) - t210;
	t224 = sin(t204);
	t223 = pkin(2) * sin(qJ(2));
	t211 = sin(pkin(10));
	t212 = sin(pkin(6));
	t221 = t211 * t212;
	t213 = cos(pkin(10));
	t220 = t212 * t213;
	t214 = cos(pkin(6));
	t219 = t226 * t212 - t214 * t223;
	t215 = -qJ(5) - pkin(8);
	t209 = qJ(4) + pkin(12);
	t208 = cos(t210);
	t207 = cos(t209);
	t206 = sin(t210);
	t205 = sin(t209);
	t202 = cos(qJ(2)) * pkin(2) + pkin(1);
	t201 = cos(qJ(4)) * pkin(4) + pkin(3);
	t200 = cos(t203);
	t198 = cos(t204) / 0.2e1;
	t195 = t198 - t200 / 0.2e1;
	t194 = t200 / 0.2e1 + t198;
	t193 = t224 / 0.2e1 + t225;
	t192 = t225 - t224 / 0.2e1;
	t191 = -t211 * t192 + t213 * t208;
	t190 = t211 * t194 + t213 * t206;
	t189 = t213 * t192 + t211 * t208;
	t188 = -t213 * t194 + t211 * t206;
	t1 = [t191 * t207 + t205 * t221, -t191 * t205 + t207 * t221, t190, -t190 * t215 + t191 * t201 + t213 * t202 + t219 * t211 + 0; t189 * t207 - t205 * t220, -t189 * t205 - t207 * t220, t188, -t188 * t215 + t189 * t201 + t211 * t202 - t219 * t213 + 0; t195 * t207 + t214 * t205, -t195 * t205 + t214 * t207, -t193, t193 * t215 + t195 * t201 + t212 * t223 + t226 * t214 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:58:39
	% EndTime: 2020-11-04 20:58:39
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (165->62), mult. (226->102), div. (0->0), fcn. (268->18), ass. (0->46)
	t247 = sin(pkin(12));
	t251 = cos(pkin(12));
	t236 = t251 * pkin(5) + t247 * pkin(9) + pkin(4);
	t237 = -t247 * pkin(5) + t251 * pkin(9);
	t258 = sin(qJ(4));
	t261 = cos(qJ(4));
	t264 = t236 * t261 + t237 * t258;
	t245 = qJ(4) + pkin(12);
	t242 = cos(t245);
	t246 = qJ(2) + pkin(11);
	t243 = cos(t246);
	t278 = t242 * t243;
	t249 = sin(pkin(10));
	t254 = cos(pkin(6));
	t277 = t249 * t254;
	t250 = sin(pkin(6));
	t276 = t250 * t242;
	t253 = cos(pkin(10));
	t275 = t250 * t253;
	t257 = sin(qJ(6));
	t274 = t250 * t257;
	t260 = cos(qJ(6));
	t273 = t250 * t260;
	t252 = cos(pkin(11));
	t272 = t252 * t249;
	t271 = t253 * t254;
	t270 = t254 * t257;
	t269 = t254 * t260;
	t248 = sin(pkin(11));
	t255 = qJ(5) + pkin(8);
	t234 = t252 * pkin(3) + t255 * t248 + pkin(2);
	t268 = t242 * t270;
	t267 = t242 * t269;
	t266 = t257 * t278;
	t265 = t260 * t278;
	t263 = t236 * t258 - t237 * t261 + pkin(7) + qJ(3);
	t262 = cos(qJ(2));
	t259 = sin(qJ(2));
	t241 = sin(t246);
	t240 = sin(t245);
	t239 = t255 * t252;
	t235 = -t248 * pkin(3) + t239;
	t229 = t254 * t240 + t241 * t276;
	t228 = t240 * t273 + t243 * t270;
	t227 = t240 * t274 - t243 * t269;
	t1 = [(-t249 * t267 + t257 * t253) * t241 + t253 * t265 + t249 * t228, (t249 * t268 + t260 * t253) * t241 - t253 * t266 - t249 * t227, (-t241 * t277 + t253 * t243) * t240 - t249 * t276, (t253 * t234 - t264 * (t248 * t277 - t253 * t252)) * t262 + (t253 * t235 - t264 * (t253 * t248 + t254 * t272)) * t259 + t253 * pkin(1) + 0 + ((-t234 * t259 + t235 * t262) * t254 + t263 * t250) * t249; (t257 * t249 + t253 * t267) * t241 + t249 * t265 - t253 * t228, (t260 * t249 - t253 * t268) * t241 - t249 * t266 + t253 * t227, (t241 * t271 + t249 * t243) * t240 + t242 * t275, (-t235 * t271 + t249 * t234 + t264 * (t248 * t271 + t272)) * t262 + (t234 * t271 + t249 * t235 + t264 * (-t249 * t248 + t252 * t271)) * t259 + t249 * pkin(1) + 0 - t263 * t275; t229 * t260 - t243 * t274, -t229 * t257 - t243 * t273, t250 * t241 * t240 - t254 * t242, qJ(1) + 0 + ((t264 * t252 + t234) * t259 + (-t239 - (-pkin(3) - t264) * t248) * t262) * t250 + t263 * t254; 0, 0, 0, 1;];
	Tc_mdh = t1;
end