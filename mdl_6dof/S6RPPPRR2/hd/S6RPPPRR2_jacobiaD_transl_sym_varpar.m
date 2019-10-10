% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPPPRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2,theta4]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:27
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPPPRR2_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR2_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR2_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPPRR2_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPPRR2_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR2_jacobiaD_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:27:33
	% EndTime: 2019-10-09 23:27:33
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:27:33
	% EndTime: 2019-10-09 23:27:33
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
	% StartTime: 2019-10-09 23:27:33
	% EndTime: 2019-10-09 23:27:33
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
	% StartTime: 2019-10-09 23:27:33
	% EndTime: 2019-10-09 23:27:33
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
	% StartTime: 2019-10-09 23:27:33
	% EndTime: 2019-10-09 23:27:33
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (37->13), mult. (40->16), div. (0->0), fcn. (26->6), ass. (0->8)
	t16 = qJ(1) + pkin(9);
	t14 = sin(t16);
	t21 = qJD(1) * t14;
	t20 = -pkin(2) - r_i_i_C(3) - qJ(4);
	t19 = r_i_i_C(1) * sin(pkin(10)) + r_i_i_C(2) * cos(pkin(10)) + qJ(3);
	t15 = cos(t16);
	t13 = qJD(1) * t15;
	t1 = [t15 * qJD(3) - t14 * qJD(4) + (-cos(qJ(1)) * pkin(1) + t20 * t15 - t19 * t14) * qJD(1), 0, t13, -t21, 0, 0; t14 * qJD(3) + t15 * qJD(4) + (-sin(qJ(1)) * pkin(1) + t19 * t15 + t20 * t14) * qJD(1), 0, t21, t13, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:27:33
	% EndTime: 2019-10-09 23:27:34
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (84->24), mult. (86->34), div. (0->0), fcn. (56->7), ass. (0->15)
	t30 = pkin(10) + qJ(5);
	t26 = sin(t30);
	t28 = cos(t30);
	t41 = (r_i_i_C(1) * t28 - r_i_i_C(2) * t26) * qJD(5);
	t31 = qJ(1) + pkin(9);
	t27 = sin(t31);
	t40 = qJD(1) * t27;
	t29 = cos(t31);
	t25 = qJD(1) * t29;
	t39 = qJD(5) * t27;
	t38 = qJD(5) * t29;
	t37 = -pkin(2) - r_i_i_C(3) - pkin(7) - qJ(4);
	t35 = pkin(4) * sin(pkin(10)) + r_i_i_C(1) * t26 + r_i_i_C(2) * t28 + qJ(3);
	t34 = qJD(3) + t41;
	t1 = [-t27 * qJD(4) + t34 * t29 + (-cos(qJ(1)) * pkin(1) + t37 * t29 - t35 * t27) * qJD(1), 0, t25, -t40, (-t26 * t25 - t28 * t39) * r_i_i_C(2) + (t28 * t25 - t26 * t39) * r_i_i_C(1), 0; t29 * qJD(4) + t34 * t27 + (-sin(qJ(1)) * pkin(1) + t37 * t27 + t35 * t29) * qJD(1), 0, t40, t25, (-t26 * t40 + t28 * t38) * r_i_i_C(2) + (t26 * t38 + t28 * t40) * r_i_i_C(1), 0; 0, 0, 0, 0, -t41, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:27:35
	% EndTime: 2019-10-09 23:27:35
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (259->47), mult. (292->74), div. (0->0), fcn. (227->9), ass. (0->34)
	t217 = pkin(10) + qJ(5);
	t213 = sin(t217);
	t215 = cos(t217);
	t239 = pkin(8) + r_i_i_C(3);
	t232 = t239 * t215;
	t244 = t232 - pkin(4) * sin(pkin(10)) - pkin(5) * t213 - qJ(3);
	t221 = sin(qJ(6));
	t222 = cos(qJ(6));
	t228 = r_i_i_C(1) * t222 - r_i_i_C(2) * t221 + pkin(5);
	t233 = t239 * t213;
	t243 = t228 * t215 + t233;
	t230 = qJD(1) * t213 + qJD(6);
	t242 = t221 * t230;
	t241 = t222 * t230;
	t238 = -pkin(2) - pkin(7) - qJ(4);
	t218 = qJ(1) + pkin(9);
	t214 = sin(t218);
	t237 = qJD(1) * t214;
	t216 = cos(t218);
	t212 = qJD(1) * t216;
	t236 = qJD(5) * t221;
	t235 = qJD(5) * t222;
	t234 = qJD(6) * t215;
	t231 = -qJD(6) * t213 - qJD(1);
	t229 = r_i_i_C(1) * t221 + r_i_i_C(2) * t222;
	t226 = t229 * qJD(6);
	t225 = qJD(3) + (pkin(5) * t215 + t233) * qJD(5);
	t224 = -t215 * t236 + t231 * t222;
	t223 = t215 * t235 + t231 * t221;
	t211 = t223 * t214 + t216 * t241;
	t210 = t224 * t214 - t216 * t242;
	t209 = -t214 * t241 + t223 * t216;
	t208 = t214 * t242 + t224 * t216;
	t1 = [t209 * r_i_i_C(1) + t208 * r_i_i_C(2) - t214 * qJD(4) + t225 * t216 + (-cos(qJ(1)) * pkin(1) + t238 * t216 + t244 * t214) * qJD(1), 0, t212, -t237, t243 * t212 + (-t229 * t234 + (-t228 * t213 + t232) * qJD(5)) * t214, t210 * r_i_i_C(1) - t211 * r_i_i_C(2); t211 * r_i_i_C(1) + t210 * r_i_i_C(2) + t216 * qJD(4) + t225 * t214 + (-sin(qJ(1)) * pkin(1) + t238 * t214 - t244 * t216) * qJD(1), 0, t237, t212, (t228 * t216 * qJD(5) + t239 * t237) * t213 + (t228 * t237 + (-t239 * qJD(5) + t226) * t216) * t215, -t208 * r_i_i_C(1) + t209 * r_i_i_C(2); 0, 0, 0, 0, -t243 * qJD(5) + t213 * t226, (t213 * t235 + t221 * t234) * r_i_i_C(2) + (t213 * t236 - t222 * t234) * r_i_i_C(1);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end