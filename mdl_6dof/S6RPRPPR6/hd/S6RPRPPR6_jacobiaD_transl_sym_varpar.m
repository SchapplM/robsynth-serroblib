% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPRPPR6
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4,theta5]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:24
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRPPR6_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR6_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR6_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPPR6_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPPR6_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR6_jacobiaD_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:24:00
	% EndTime: 2019-10-10 00:24:00
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:24:00
	% EndTime: 2019-10-10 00:24:00
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
	% StartTime: 2019-10-10 00:24:00
	% EndTime: 2019-10-10 00:24:00
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
	% StartTime: 2019-10-10 00:24:00
	% EndTime: 2019-10-10 00:24:01
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
	% StartTime: 2019-10-10 00:24:00
	% EndTime: 2019-10-10 00:24:01
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
	% StartTime: 2019-10-10 00:24:02
	% EndTime: 2019-10-10 00:24:02
	% DurationCPUTime: 0.21s
	% Computational Cost: add. (128->25), mult. (216->33), div. (0->0), fcn. (160->8), ass. (0->20)
	t174 = qJ(3) + pkin(9);
	t171 = sin(t174);
	t172 = cos(t174);
	t175 = sin(pkin(10));
	t176 = cos(pkin(10));
	t187 = r_i_i_C(1) * t176 - r_i_i_C(2) * t175 + pkin(4);
	t192 = r_i_i_C(3) + qJ(5);
	t185 = t187 * t171 - t192 * t172 + sin(qJ(3)) * pkin(3);
	t199 = qJD(4) + (qJ(2) + t185) * qJD(1);
	t197 = qJD(3) * t185 - qJD(5) * t171;
	t184 = t192 * t171 + t187 * t172 + cos(qJ(3)) * pkin(3);
	t196 = -t184 * qJD(3) + t172 * qJD(5);
	t179 = sin(qJ(1));
	t191 = qJD(1) * t179;
	t181 = cos(qJ(1));
	t173 = qJD(1) * t181;
	t190 = qJD(3) * t171;
	t183 = qJD(1) * t184;
	t182 = qJD(2) + (-r_i_i_C(1) * t175 - r_i_i_C(2) * t176 - pkin(1) - pkin(7) - qJ(4)) * qJD(1) - t196;
	t1 = [-t179 * t199 + t182 * t181, t173, -t179 * t197 + t181 * t183, -t191, -t172 * t173 + t179 * t190, 0; t182 * t179 + t181 * t199, t191, t179 * t183 + t181 * t197, t173, -t172 * t191 - t181 * t190, 0; 0, 0, t196, 0, qJD(3) * t172, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:24:02
	% EndTime: 2019-10-10 00:24:02
	% DurationCPUTime: 0.29s
	% Computational Cost: add. (269->50), mult. (341->75), div. (0->0), fcn. (269->10), ass. (0->38)
	t221 = qJ(3) + pkin(9);
	t216 = sin(t221);
	t218 = cos(t221);
	t251 = r_i_i_C(3) + pkin(8) + qJ(5);
	t234 = t251 * t218 - sin(qJ(3)) * pkin(3);
	t214 = cos(pkin(10)) * pkin(5) + pkin(4);
	t220 = pkin(10) + qJ(6);
	t215 = sin(t220);
	t217 = cos(t220);
	t235 = r_i_i_C(1) * t217 - r_i_i_C(2) * t215 + t214;
	t263 = -(-t235 * t216 + t234) * qJD(3) - qJD(5) * t216;
	t262 = -qJD(4) + (-t214 * t216 - qJ(2) + t234) * qJD(1);
	t232 = t251 * t216 + cos(qJ(3)) * pkin(3);
	t259 = t235 * t218 + t232;
	t226 = sin(qJ(1));
	t239 = qJD(1) * t216 + qJD(6);
	t228 = cos(qJ(1));
	t246 = qJD(3) * t228;
	t255 = -t218 * t246 + t239 * t226;
	t247 = qJD(3) * t226;
	t254 = t218 * t247 + t239 * t228;
	t249 = qJD(1) * t226;
	t219 = qJD(1) * t228;
	t248 = qJD(3) * t216;
	t244 = qJD(6) * t218;
	t243 = t218 * qJD(5);
	t240 = -qJD(6) * t216 - qJD(1);
	t238 = r_i_i_C(1) * t215 + r_i_i_C(2) * t217;
	t237 = t240 * t226;
	t236 = t240 * t228;
	t231 = qJD(6) * t238;
	t230 = qJD(1) * t259;
	t229 = -t243 + qJD(2) + (-pkin(5) * sin(pkin(10)) - pkin(1) - qJ(4) - pkin(7)) * qJD(1) + (t214 * t218 + t232) * qJD(3);
	t213 = t215 * t237 + t254 * t217;
	t212 = -t254 * t215 + t217 * t237;
	t211 = t215 * t236 - t255 * t217;
	t210 = t255 * t215 + t217 * t236;
	t1 = [t211 * r_i_i_C(1) + t210 * r_i_i_C(2) + t262 * t226 + t229 * t228, t219, t228 * t230 + (-t238 * t244 - t263) * t226, -t249, t216 * t247 - t218 * t219, t212 * r_i_i_C(1) - t213 * r_i_i_C(2); t213 * r_i_i_C(1) + t212 * r_i_i_C(2) + t229 * t226 - t262 * t228, t249, t226 * t230 + (t218 * t231 + t263) * t228, t219, -t216 * t246 - t218 * t249, -t210 * r_i_i_C(1) + t211 * r_i_i_C(2); 0, 0, -t259 * qJD(3) + t216 * t231 + t243, 0, qJD(3) * t218, (t215 * t244 + t217 * t248) * r_i_i_C(2) + (t215 * t248 - t217 * t244) * r_i_i_C(1);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end