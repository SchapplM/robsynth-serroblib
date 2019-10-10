% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPRPRP7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:39
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRPRP7_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP7_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP7_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRP7_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRP7_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP7_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:39:23
	% EndTime: 2019-10-10 00:39:23
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:39:23
	% EndTime: 2019-10-10 00:39:23
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
	% StartTime: 2019-10-10 00:39:23
	% EndTime: 2019-10-10 00:39:23
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (8->6), mult. (20->10), div. (0->0), fcn. (12->2), ass. (0->5)
	t10 = -pkin(1) + r_i_i_C(2);
	t9 = r_i_i_C(3) + qJ(2);
	t8 = cos(qJ(1));
	t7 = sin(qJ(1));
	t1 = [t8 * qJD(2) + (t10 * t8 - t9 * t7) * qJD(1), qJD(1) * t8, 0, 0, 0, 0; t7 * qJD(2) + (t10 * t7 + t9 * t8) * qJD(1), qJD(1) * t7, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:39:23
	% EndTime: 2019-10-10 00:39:23
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (23->16), mult. (72->29), div. (0->0), fcn. (46->4), ass. (0->13)
	t17 = sin(qJ(3));
	t19 = cos(qJ(3));
	t29 = (r_i_i_C(1) * t19 - r_i_i_C(2) * t17) * qJD(3);
	t18 = sin(qJ(1));
	t28 = qJD(1) * t18;
	t20 = cos(qJ(1));
	t27 = qJD(1) * t20;
	t26 = qJD(3) * t18;
	t25 = qJD(3) * t20;
	t24 = -pkin(1) - pkin(7) - r_i_i_C(3);
	t22 = r_i_i_C(1) * t17 + r_i_i_C(2) * t19 + qJ(2);
	t21 = qJD(2) + t29;
	t1 = [t21 * t20 + (-t18 * t22 + t20 * t24) * qJD(1), t27, (-t17 * t27 - t19 * t26) * r_i_i_C(2) + (-t17 * t26 + t19 * t27) * r_i_i_C(1), 0, 0, 0; t21 * t18 + (t18 * t24 + t20 * t22) * qJD(1), t28, (-t17 * t28 + t19 * t25) * r_i_i_C(2) + (t17 * t25 + t19 * t28) * r_i_i_C(1), 0, 0, 0; 0, 0, -t29, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:39:23
	% EndTime: 2019-10-10 00:39:23
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (55->18), mult. (102->25), div. (0->0), fcn. (67->6), ass. (0->16)
	t25 = qJ(3) + pkin(9);
	t22 = sin(t25);
	t23 = cos(t25);
	t35 = sin(qJ(3)) * pkin(3) + r_i_i_C(1) * t22 + r_i_i_C(2) * t23;
	t39 = qJD(3) * t35;
	t34 = cos(qJ(3)) * pkin(3) + r_i_i_C(1) * t23 - r_i_i_C(2) * t22;
	t38 = t34 * qJD(3);
	t28 = sin(qJ(1));
	t37 = qJD(1) * t28;
	t36 = -pkin(1) - r_i_i_C(3) - qJ(4) - pkin(7);
	t33 = qJ(2) + t35;
	t32 = qJD(1) * t34;
	t31 = qJD(2) + t38;
	t30 = cos(qJ(1));
	t24 = qJD(1) * t30;
	t1 = [-t28 * qJD(4) + t31 * t30 + (-t28 * t33 + t30 * t36) * qJD(1), t24, -t28 * t39 + t30 * t32, -t37, 0, 0; t30 * qJD(4) + t31 * t28 + (t28 * t36 + t30 * t33) * qJD(1), t37, t28 * t32 + t30 * t39, t24, 0, 0; 0, 0, -t38, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:39:25
	% EndTime: 2019-10-10 00:39:25
	% DurationCPUTime: 0.27s
	% Computational Cost: add. (172->44), mult. (308->72), div. (0->0), fcn. (238->8), ass. (0->34)
	t211 = qJ(3) + pkin(9);
	t208 = sin(t211);
	t209 = cos(t211);
	t239 = pkin(8) + r_i_i_C(3);
	t224 = t239 * t209 - sin(qJ(3)) * pkin(3);
	t213 = sin(qJ(5));
	t216 = cos(qJ(5));
	t226 = r_i_i_C(1) * t216 - r_i_i_C(2) * t213 + pkin(4);
	t248 = (-t226 * t208 + t224) * qJD(3);
	t247 = -pkin(4) * t208 - qJ(2) + t224;
	t222 = t239 * t208 + cos(qJ(3)) * pkin(3);
	t245 = t226 * t209 + t222;
	t218 = cos(qJ(1));
	t228 = qJD(1) * t208 + qJD(5);
	t242 = t228 * t218;
	t215 = sin(qJ(1));
	t233 = qJD(3) * t209;
	t240 = t228 * t215 - t218 * t233;
	t235 = -pkin(1) - qJ(4) - pkin(7);
	t234 = qJD(1) * t215;
	t232 = qJD(3) * t216;
	t231 = qJD(5) * t209;
	t229 = -qJD(5) * t208 - qJD(1);
	t227 = r_i_i_C(1) * t213 + r_i_i_C(2) * t216;
	t225 = t229 * t218;
	t221 = qJD(5) * t227;
	t220 = qJD(2) + (pkin(4) * t209 + t222) * qJD(3);
	t219 = qJD(1) * t245;
	t210 = qJD(1) * t218;
	t207 = t216 * t242 + (t209 * t232 + t229 * t213) * t215;
	t206 = t229 * t216 * t215 + (-t215 * t233 - t242) * t213;
	t205 = t213 * t225 - t240 * t216;
	t204 = t240 * t213 + t216 * t225;
	t1 = [t205 * r_i_i_C(1) + t204 * r_i_i_C(2) - t215 * qJD(4) + t220 * t218 + (t247 * t215 + t235 * t218) * qJD(1), t210, t218 * t219 + (-t227 * t231 + t248) * t215, -t234, t206 * r_i_i_C(1) - t207 * r_i_i_C(2), 0; t207 * r_i_i_C(1) + t206 * r_i_i_C(2) + t218 * qJD(4) + t220 * t215 + (t235 * t215 - t247 * t218) * qJD(1), t234, t215 * t219 + (t209 * t221 - t248) * t218, t210, -t204 * r_i_i_C(1) + t205 * r_i_i_C(2), 0; 0, 0, -t245 * qJD(3) + t208 * t221, 0, (t208 * t232 + t213 * t231) * r_i_i_C(2) + (qJD(3) * t208 * t213 - t216 * t231) * r_i_i_C(1), 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:39:25
	% EndTime: 2019-10-10 00:39:25
	% DurationCPUTime: 0.33s
	% Computational Cost: add. (241->54), mult. (408->79), div. (0->0), fcn. (317->8), ass. (0->42)
	t219 = cos(qJ(5));
	t247 = t219 * pkin(5);
	t209 = pkin(4) + t247;
	t213 = qJ(3) + pkin(9);
	t210 = sin(t213);
	t211 = cos(t213);
	t246 = r_i_i_C(3) + qJ(6) + pkin(8);
	t228 = t246 * t211 - sin(qJ(3)) * pkin(3);
	t263 = (-t209 * t210 - qJ(2) + t228) * qJD(1) - qJD(5) * t247 - qJD(4);
	t216 = sin(qJ(5));
	t249 = r_i_i_C(2) * t216;
	t230 = r_i_i_C(1) * t219 + t209 - t249;
	t262 = -(-t230 * t210 + t228) * qJD(3) - qJD(6) * t210;
	t218 = sin(qJ(1));
	t233 = qJD(1) * t210 + qJD(5);
	t221 = cos(qJ(1));
	t240 = qJD(3) * t221;
	t254 = -t211 * t240 + t233 * t218;
	t261 = t254 * t216;
	t226 = t246 * t210 + cos(qJ(3)) * pkin(3);
	t258 = t230 * t211 + t226;
	t253 = pkin(5) + r_i_i_C(1);
	t229 = r_i_i_C(2) * t219 + t253 * t216;
	t250 = pkin(5) * t216;
	t244 = t219 * t221;
	t243 = qJD(1) * t218;
	t212 = qJD(1) * t221;
	t242 = qJD(3) * t211;
	t241 = qJD(3) * t218;
	t239 = qJD(5) * t210;
	t238 = qJD(5) * t211;
	t236 = t211 * qJD(6);
	t234 = qJD(1) + t239;
	t231 = t234 * t221;
	t225 = qJD(5) * t229;
	t223 = qJD(1) * t258;
	t207 = -t234 * t219 * t218 + (-t211 * t241 - t233 * t221) * t216;
	t222 = -t239 * t250 - t236 + qJD(2) + (-pkin(1) - qJ(4) - pkin(7) - t250) * qJD(1) + (t209 * t211 + t226) * qJD(3);
	t208 = t233 * t244 + (-t234 * t216 + t219 * t242) * t218;
	t206 = -t216 * t231 - t219 * t254;
	t205 = -t219 * t231 + t261;
	t1 = [t206 * r_i_i_C(1) + t205 * r_i_i_C(2) + t263 * t218 + t222 * t221, t212, t221 * t223 + (-t229 * t238 - t262) * t218, -t243, -t208 * r_i_i_C(2) + t253 * t207, t210 * t241 - t211 * t212; t208 * r_i_i_C(1) + t207 * r_i_i_C(2) + t222 * t218 - t263 * t221, t243, t218 * t223 + (t211 * t225 + t262) * t221, t212, -t205 * r_i_i_C(1) + t206 * r_i_i_C(2) + (t234 * t244 - t261) * pkin(5), -t210 * t240 - t211 * t243; 0, 0, -t258 * qJD(3) + t210 * t225 + t236, 0, (-t253 * t219 + t249) * t238 + t229 * t210 * qJD(3), t242;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end