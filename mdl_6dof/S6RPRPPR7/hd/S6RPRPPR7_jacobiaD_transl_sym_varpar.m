% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPRPPR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:25
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRPPR7_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR7_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR7_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPPR7_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPPR7_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR7_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:25:42
	% EndTime: 2019-10-10 00:25:42
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:25:42
	% EndTime: 2019-10-10 00:25:42
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
	% StartTime: 2019-10-10 00:25:42
	% EndTime: 2019-10-10 00:25:42
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
	% StartTime: 2019-10-10 00:25:42
	% EndTime: 2019-10-10 00:25:42
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
	% StartTime: 2019-10-10 00:25:42
	% EndTime: 2019-10-10 00:25:42
	% DurationCPUTime: 0.11s
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
	t1 = [-t28 * qJD(4) + t31 * t30 + (-t33 * t28 + t36 * t30) * qJD(1), t24, -t28 * t39 + t30 * t32, -t37, 0, 0; t30 * qJD(4) + t31 * t28 + (t36 * t28 + t33 * t30) * qJD(1), t37, t28 * t32 + t30 * t39, t24, 0, 0; 0, 0, -t38, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:25:43
	% EndTime: 2019-10-10 00:25:43
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (108->23), mult. (172->29), div. (0->0), fcn. (119->6), ass. (0->18)
	t148 = qJ(3) + pkin(9);
	t145 = sin(t148);
	t146 = cos(t148);
	t163 = r_i_i_C(3) + qJ(5);
	t166 = pkin(4) - r_i_i_C(2);
	t157 = t166 * t145 - t163 * t146 + sin(qJ(3)) * pkin(3);
	t171 = qJD(4) + (qJ(2) + t157) * qJD(1);
	t169 = t157 * qJD(3) - qJD(5) * t145;
	t156 = t163 * t145 + t166 * t146 + cos(qJ(3)) * pkin(3);
	t168 = -t156 * qJD(3) + t146 * qJD(5);
	t151 = sin(qJ(1));
	t162 = qJD(1) * t151;
	t153 = cos(qJ(1));
	t147 = qJD(1) * t153;
	t161 = qJD(3) * t145;
	t155 = qJD(1) * t156;
	t154 = qJD(2) + (-pkin(1) - r_i_i_C(1) - qJ(4) - pkin(7)) * qJD(1) - t168;
	t1 = [-t171 * t151 + t154 * t153, t147, -t169 * t151 + t153 * t155, -t162, -t146 * t147 + t151 * t161, 0; t154 * t151 + t171 * t153, t162, t151 * t155 + t169 * t153, t147, -t146 * t162 - t153 * t161, 0; 0, 0, t168, 0, qJD(3) * t146, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:25:43
	% EndTime: 2019-10-10 00:25:44
	% DurationCPUTime: 0.32s
	% Computational Cost: add. (209->47), mult. (356->70), div. (0->0), fcn. (275->8), ass. (0->35)
	t211 = qJ(3) + pkin(9);
	t208 = sin(t211);
	t209 = cos(t211);
	t213 = sin(qJ(6));
	t216 = cos(qJ(6));
	t222 = qJD(5) + (r_i_i_C(1) * t216 - r_i_i_C(2) * t213) * qJD(6);
	t233 = pkin(4) + pkin(8) + r_i_i_C(3);
	t225 = t233 * t208 + sin(qJ(3)) * pkin(3);
	t226 = r_i_i_C(1) * t213 + r_i_i_C(2) * t216 + qJ(5);
	t250 = -t222 * t208 + (-t226 * t209 + t225) * qJD(3);
	t249 = qJD(4) + (-qJ(5) * t209 + qJ(2) + t225) * qJD(1);
	t224 = t233 * t209 + cos(qJ(3)) * pkin(3);
	t243 = t226 * t208 + t224;
	t218 = cos(qJ(1));
	t239 = t216 * t218;
	t215 = sin(qJ(1));
	t238 = qJD(1) * t215;
	t210 = qJD(1) * t218;
	t237 = qJD(3) * t208;
	t236 = qJD(3) * t209;
	t235 = qJD(3) * t216;
	t234 = qJD(6) * t208;
	t232 = t215 * t237;
	t231 = t218 * t237;
	t230 = qJD(6) * t209 + qJD(1);
	t229 = qJD(1) * t209 + qJD(6);
	t227 = t230 * t213;
	t221 = t229 * t215 + t231;
	t220 = qJD(1) * t243;
	t219 = -t209 * qJD(5) + qJD(2) + (-pkin(1) - pkin(5) - qJ(4) - pkin(7)) * qJD(1) + (qJ(5) * t208 + t224) * qJD(3);
	t207 = t221 * t213 - t230 * t239;
	t206 = t221 * t216 + t218 * t227;
	t205 = t230 * t216 * t215 + (t229 * t218 - t232) * t213;
	t204 = -t229 * t239 + (t208 * t235 + t227) * t215;
	t1 = [t207 * r_i_i_C(1) + t206 * r_i_i_C(2) - t249 * t215 + t219 * t218, t210, -t250 * t215 + t218 * t220, -t238, -t209 * t210 + t232, t204 * r_i_i_C(1) + t205 * r_i_i_C(2); -t205 * r_i_i_C(1) + t204 * r_i_i_C(2) + t219 * t215 + t249 * t218, t238, t215 * t220 + t250 * t218, t210, -t209 * t238 - t231, -t206 * r_i_i_C(1) + t207 * r_i_i_C(2); 0, 0, -t243 * qJD(3) + t222 * t209, 0, t236, (-t213 * t236 - t216 * t234) * r_i_i_C(2) + (t209 * t235 - t213 * t234) * r_i_i_C(1);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end