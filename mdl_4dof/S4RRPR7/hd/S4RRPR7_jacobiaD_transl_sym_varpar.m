% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S4RRPR7
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% JaD_transl [3x4]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 13:41
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S4RRPR7_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),uint8(0),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR7_jacobiaD_transl_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR7_jacobiaD_transl_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S4RRPR7_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4RRPR7_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR7_jacobiaD_transl_sym_varpar: pkin has to be [7x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 13:41:28
	% EndTime: 2019-12-29 13:41:28
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 13:41:26
	% EndTime: 2019-12-29 13:41:26
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0; 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 13:41:28
	% EndTime: 2019-12-29 13:41:28
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (19->15), mult. (64->29), div. (0->0), fcn. (40->4), ass. (0->13)
	t28 = pkin(5) + r_i_i_C(3);
	t18 = sin(qJ(1));
	t27 = qJD(1) * t18;
	t20 = cos(qJ(1));
	t26 = qJD(1) * t20;
	t25 = qJD(2) * t18;
	t24 = qJD(2) * t20;
	t17 = sin(qJ(2));
	t19 = cos(qJ(2));
	t23 = r_i_i_C(1) * t17 + r_i_i_C(2) * t19;
	t22 = -r_i_i_C(1) * t19 + r_i_i_C(2) * t17 - pkin(1);
	t21 = t23 * qJD(2);
	t1 = [t23 * t25 + (-t18 * t28 + t22 * t20) * qJD(1), (t17 * t24 + t19 * t27) * r_i_i_C(2) + (t17 * t27 - t19 * t24) * r_i_i_C(1), 0, 0; -t20 * t21 + (t22 * t18 + t20 * t28) * qJD(1), (t17 * t25 - t19 * t26) * r_i_i_C(2) + (-t17 * t26 - t19 * t25) * r_i_i_C(1), 0, 0; 0, -t21, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 13:41:26
	% EndTime: 2019-12-29 13:41:27
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (50->16), mult. (90->25), div. (0->0), fcn. (59->6), ass. (0->14)
	t27 = qJ(2) + pkin(7);
	t25 = sin(t27);
	t26 = cos(t27);
	t42 = -r_i_i_C(1) * t26 + r_i_i_C(2) * t25 - cos(qJ(2)) * pkin(2);
	t40 = r_i_i_C(3) + qJ(3) + pkin(5);
	t32 = cos(qJ(1));
	t39 = qJD(1) * t32;
	t37 = -pkin(1) + t42;
	t36 = sin(qJ(2)) * pkin(2) + r_i_i_C(1) * t25 + r_i_i_C(2) * t26;
	t30 = sin(qJ(1));
	t35 = t36 * t30;
	t34 = qJD(2) * t42;
	t33 = t36 * qJD(2);
	t1 = [t32 * qJD(3) + qJD(2) * t35 + (-t40 * t30 + t37 * t32) * qJD(1), qJD(1) * t35 + t32 * t34, t39, 0; t30 * qJD(3) - t32 * t33 + (t37 * t30 + t40 * t32) * qJD(1), t30 * t34 - t36 * t39, qJD(1) * t30, 0; 0, -t33, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 13:41:30
	% EndTime: 2019-12-29 13:41:31
	% DurationCPUTime: 0.32s
	% Computational Cost: add. (167->41), mult. (296->69), div. (0->0), fcn. (230->8), ass. (0->34)
	t212 = qJ(2) + pkin(7);
	t210 = sin(t212);
	t211 = cos(t212);
	t243 = pkin(6) + r_i_i_C(3);
	t225 = t243 * t211 - sin(qJ(2)) * pkin(2);
	t246 = -pkin(3) * t210 + t225;
	t214 = sin(qJ(4));
	t217 = cos(qJ(4));
	t227 = r_i_i_C(1) * t217 - r_i_i_C(2) * t214 + pkin(3);
	t221 = -t227 * t210 + t225;
	t244 = -t243 * t210 - cos(qJ(2)) * pkin(2);
	t219 = cos(qJ(1));
	t239 = t217 * t219;
	t216 = sin(qJ(1));
	t238 = qJD(1) * t216;
	t237 = qJD(1) * t219;
	t236 = qJD(2) * t216;
	t235 = qJD(2) * t217;
	t234 = qJD(2) * t219;
	t233 = qJD(4) * t210;
	t232 = qJD(4) * t211;
	t230 = -qJD(1) + t232;
	t229 = qJD(1) * t211 - qJD(4);
	t228 = r_i_i_C(1) * t214 + r_i_i_C(2) * t217;
	t226 = t230 * t214;
	t223 = -pkin(3) * t211 - pkin(1) + t244;
	t222 = t210 * t234 + t229 * t216;
	t220 = t228 * t233 + (-t227 * t211 + t244) * qJD(2);
	t213 = -qJ(3) - pkin(5);
	t208 = -t229 * t239 + (t210 * t235 + t226) * t216;
	t207 = t230 * t217 * t216 + (-t210 * t236 + t229 * t219) * t214;
	t206 = t222 * t217 + t219 * t226;
	t205 = t222 * t214 - t230 * t239;
	t1 = [t208 * r_i_i_C(1) + t207 * r_i_i_C(2) + t219 * qJD(3) - t246 * t236 + (t213 * t216 + t223 * t219) * qJD(1), t220 * t219 - t221 * t238, t237, t205 * r_i_i_C(1) + t206 * r_i_i_C(2); -t206 * r_i_i_C(1) + t205 * r_i_i_C(2) + t216 * qJD(3) + t246 * t234 + (-t213 * t219 + t223 * t216) * qJD(1), t220 * t216 + t221 * t237, t238, -t207 * r_i_i_C(1) + t208 * r_i_i_C(2); 0, t221 * qJD(2) - t228 * t232, 0, (-t211 * t235 + t214 * t233) * r_i_i_C(2) + (-qJD(2) * t211 * t214 - t217 * t233) * r_i_i_C(1);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,4);
end