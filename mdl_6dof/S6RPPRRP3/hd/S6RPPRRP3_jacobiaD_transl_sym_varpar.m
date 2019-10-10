% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPPRRP3
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:51
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPPRRP3_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP3_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP3_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRRP3_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRRP3_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP3_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:51:14
	% EndTime: 2019-10-09 23:51:14
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:51:14
	% EndTime: 2019-10-09 23:51:14
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:51:14
	% EndTime: 2019-10-09 23:51:14
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->5), mult. (12->8), div. (0->0), fcn. (6->4), ass. (0->4)
	t32 = qJ(1) + pkin(9);
	t31 = cos(t32);
	t30 = sin(t32);
	t1 = [(-cos(qJ(1)) * pkin(1) - r_i_i_C(1) * t31 + r_i_i_C(2) * t30) * qJD(1), 0, 0, 0, 0, 0; (-sin(qJ(1)) * pkin(1) - r_i_i_C(1) * t30 - r_i_i_C(2) * t31) * qJD(1), 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:51:14
	% EndTime: 2019-10-09 23:51:14
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (22->9), mult. (24->12), div. (0->0), fcn. (14->4), ass. (0->6)
	t13 = -pkin(2) + r_i_i_C(2);
	t12 = r_i_i_C(3) + qJ(3);
	t11 = qJ(1) + pkin(9);
	t10 = cos(t11);
	t9 = sin(t11);
	t1 = [t10 * qJD(3) + (-cos(qJ(1)) * pkin(1) - t12 * t9 + t13 * t10) * qJD(1), 0, qJD(1) * t10, 0, 0, 0; t9 * qJD(3) + (-sin(qJ(1)) * pkin(1) + t13 * t9 + t12 * t10) * qJD(1), 0, qJD(1) * t9, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:51:14
	% EndTime: 2019-10-09 23:51:14
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (53->19), mult. (76->33), div. (0->0), fcn. (48->6), ass. (0->14)
	t24 = sin(qJ(4));
	t25 = cos(qJ(4));
	t34 = (r_i_i_C(1) * t25 - r_i_i_C(2) * t24) * qJD(4);
	t33 = qJD(1) * t24;
	t32 = qJD(1) * t25;
	t31 = qJD(4) * t24;
	t30 = qJD(4) * t25;
	t29 = -pkin(2) - pkin(7) - r_i_i_C(3);
	t27 = r_i_i_C(1) * t24 + r_i_i_C(2) * t25 + qJ(3);
	t26 = qJD(3) + t34;
	t23 = qJ(1) + pkin(9);
	t22 = cos(t23);
	t21 = sin(t23);
	t1 = [t26 * t22 + (-cos(qJ(1)) * pkin(1) + t29 * t22 - t27 * t21) * qJD(1), 0, qJD(1) * t22, (-t21 * t30 - t22 * t33) * r_i_i_C(2) + (-t21 * t31 + t22 * t32) * r_i_i_C(1), 0, 0; t26 * t21 + (-sin(qJ(1)) * pkin(1) + t27 * t22 + t29 * t21) * qJD(1), 0, qJD(1) * t21, (-t21 * t33 + t22 * t30) * r_i_i_C(2) + (t21 * t32 + t22 * t31) * r_i_i_C(1), 0, 0; 0, 0, 0, -t34, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:51:15
	% EndTime: 2019-10-09 23:51:15
	% DurationCPUTime: 0.28s
	% Computational Cost: add. (175->40), mult. (282->67), div. (0->0), fcn. (219->8), ass. (0->33)
	t202 = sin(qJ(4));
	t201 = sin(qJ(5));
	t203 = cos(qJ(5));
	t210 = r_i_i_C(1) * t203 - r_i_i_C(2) * t201 + pkin(4);
	t204 = cos(qJ(4));
	t221 = pkin(8) + r_i_i_C(3);
	t223 = t221 * t204;
	t231 = (-t210 * t202 + t223) * qJD(4);
	t216 = t221 * t202;
	t230 = t210 * t204 + t216;
	t229 = -pkin(4) * t202 - qJ(3) + t223;
	t212 = qJD(1) * t202 + qJD(5);
	t227 = t201 * t212;
	t226 = t203 * t212;
	t222 = -pkin(2) - pkin(7);
	t220 = qJD(4) * t202;
	t219 = qJD(4) * t204;
	t218 = qJD(5) * t204;
	t213 = -qJD(5) * t202 - qJD(1);
	t211 = r_i_i_C(1) * t201 + r_i_i_C(2) * t203;
	t209 = qJD(5) * t211;
	t208 = qJD(3) + (pkin(4) * t204 + t216) * qJD(4);
	t207 = -t201 * t219 + t213 * t203;
	t206 = t213 * t201 + t203 * t219;
	t205 = qJD(1) * t230;
	t200 = qJ(1) + pkin(9);
	t199 = cos(t200);
	t198 = sin(t200);
	t197 = t206 * t198 + t199 * t226;
	t196 = t207 * t198 - t199 * t227;
	t195 = -t198 * t226 + t206 * t199;
	t194 = t198 * t227 + t207 * t199;
	t1 = [t195 * r_i_i_C(1) + t194 * r_i_i_C(2) + t208 * t199 + (-cos(qJ(1)) * pkin(1) + t222 * t199 + t229 * t198) * qJD(1), 0, qJD(1) * t199, t199 * t205 + (-t211 * t218 + t231) * t198, t196 * r_i_i_C(1) - t197 * r_i_i_C(2), 0; t197 * r_i_i_C(1) + t196 * r_i_i_C(2) + t208 * t198 + (-sin(qJ(1)) * pkin(1) + t222 * t198 - t229 * t199) * qJD(1), 0, qJD(1) * t198, t198 * t205 + (t204 * t209 - t231) * t199, -t194 * r_i_i_C(1) + t195 * r_i_i_C(2), 0; 0, 0, 0, -t230 * qJD(4) + t202 * t209, (t201 * t218 + t203 * t220) * r_i_i_C(2) + (t201 * t220 - t203 * t218) * r_i_i_C(1), 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:51:16
	% EndTime: 2019-10-09 23:51:16
	% DurationCPUTime: 0.40s
	% Computational Cost: add. (328->55), mult. (522->86), div. (0->0), fcn. (435->8), ass. (0->40)
	t237 = sin(qJ(4));
	t239 = cos(qJ(4));
	t236 = sin(qJ(5));
	t238 = cos(qJ(5));
	t264 = r_i_i_C(3) + qJ(6);
	t265 = -r_i_i_C(1) - pkin(5);
	t243 = t264 * t236 - t265 * t238 + pkin(4);
	t266 = pkin(8) + r_i_i_C(2);
	t269 = t266 * t239;
	t257 = qJD(6) * t236;
	t270 = t257 + (t265 * t236 + t264 * t238) * qJD(5);
	t274 = (t243 * t237 - t269) * qJD(4) - t270 * t239;
	t254 = t266 * t237;
	t272 = t243 * t239 + t254;
	t271 = -pkin(4) * t237 - qJ(3) + t269;
	t267 = -pkin(2) - pkin(7);
	t235 = qJ(1) + pkin(9);
	t234 = cos(t235);
	t263 = t234 * t236;
	t262 = t237 * t238;
	t261 = qJD(1) * t234;
	t260 = qJD(4) * t237;
	t259 = qJD(4) * t239;
	t258 = qJD(5) * t239;
	t256 = qJD(6) * t238;
	t252 = t234 * t262;
	t251 = t236 * t259;
	t250 = t238 * t259;
	t248 = qJD(5) * t237 + qJD(1);
	t247 = qJD(1) * t237 + qJD(5);
	t246 = t247 * t236;
	t233 = sin(t235);
	t245 = t233 * t262 + t263;
	t241 = t237 * t257 + qJD(3) + (pkin(4) * t239 + t254) * qJD(4);
	t240 = qJD(1) * t272;
	t228 = t247 * t238 * t234 + (-t248 * t236 + t250) * t233;
	t227 = t234 * t246 + (t248 * t238 + t251) * t233;
	t226 = -t234 * t250 + (t233 * t238 + t237 * t263) * qJD(5) + t245 * qJD(1);
	t225 = -qJD(5) * t252 + t233 * t246 - t234 * t251 - t238 * t261;
	t1 = [t233 * t256 + t265 * t226 - t264 * t225 + t241 * t234 + (-cos(qJ(1)) * pkin(1) + t267 * t234 + t271 * t233) * qJD(1), 0, t261, -t274 * t233 + t234 * t240, t245 * qJD(6) + t265 * t227 + t264 * t228, t227; -t234 * t256 - t265 * t228 + t264 * t227 + t241 * t233 + (-sin(qJ(1)) * pkin(1) + t267 * t233 - t271 * t234) * qJD(1), 0, qJD(1) * t233, t233 * t240 + t274 * t234, -(-t233 * t236 + t252) * qJD(6) + t264 * t226 + t265 * t225, t225; 0, 0, 0, -t272 * qJD(4) - t237 * t270, (-t264 * t258 - t265 * t260) * t236 + (-t264 * t260 + (t265 * qJD(5) + qJD(6)) * t239) * t238, -t236 * t260 + t238 * t258;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end