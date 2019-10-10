% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPPRPR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:44
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPPRPR7_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR7_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR7_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRPR7_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRPR7_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR7_jacobiaD_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:44:24
	% EndTime: 2019-10-09 23:44:24
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:44:24
	% EndTime: 2019-10-09 23:44:24
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
	% StartTime: 2019-10-09 23:44:24
	% EndTime: 2019-10-09 23:44:24
	% DurationCPUTime: 0.04s
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
	% StartTime: 2019-10-09 23:44:24
	% EndTime: 2019-10-09 23:44:24
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (15->10), mult. (36->14), div. (0->0), fcn. (24->4), ass. (0->7)
	t14 = sin(qJ(1));
	t18 = qJD(1) * t14;
	t17 = -pkin(1) - r_i_i_C(3) - qJ(3);
	t16 = r_i_i_C(1) * sin(pkin(9)) + r_i_i_C(2) * cos(pkin(9)) + qJ(2);
	t15 = cos(qJ(1));
	t11 = qJD(1) * t15;
	t1 = [t15 * qJD(2) - t14 * qJD(3) + (-t16 * t14 + t17 * t15) * qJD(1), t11, -t18, 0, 0, 0; t14 * qJD(2) + t15 * qJD(3) + (t17 * t14 + t16 * t15) * qJD(1), t18, t11, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:44:24
	% EndTime: 2019-10-09 23:44:24
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (48->21), mult. (82->32), div. (0->0), fcn. (54->5), ass. (0->14)
	t24 = pkin(9) + qJ(4);
	t21 = sin(t24);
	t22 = cos(t24);
	t36 = (r_i_i_C(1) * t22 - r_i_i_C(2) * t21) * qJD(4);
	t27 = sin(qJ(1));
	t35 = qJD(1) * t27;
	t28 = cos(qJ(1));
	t23 = qJD(1) * t28;
	t34 = qJD(4) * t27;
	t33 = qJD(4) * t28;
	t32 = -pkin(1) - r_i_i_C(3) - pkin(7) - qJ(3);
	t30 = pkin(3) * sin(pkin(9)) + r_i_i_C(1) * t21 + r_i_i_C(2) * t22 + qJ(2);
	t29 = qJD(2) + t36;
	t1 = [-t27 * qJD(3) + t29 * t28 + (-t30 * t27 + t32 * t28) * qJD(1), t23, -t35, (-t21 * t23 - t22 * t34) * r_i_i_C(2) + (-t21 * t34 + t22 * t23) * r_i_i_C(1), 0, 0; t28 * qJD(3) + t29 * t27 + (t32 * t27 + t30 * t28) * qJD(1), t35, t23, (-t21 * t35 + t22 * t33) * r_i_i_C(2) + (t21 * t33 + t22 * t35) * r_i_i_C(1), 0, 0; 0, 0, 0, -t36, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:44:25
	% EndTime: 2019-10-09 23:44:25
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (121->30), mult. (196->45), div. (0->0), fcn. (147->7), ass. (0->21)
	t171 = pkin(9) + qJ(4);
	t168 = sin(t171);
	t169 = cos(t171);
	t172 = sin(pkin(10));
	t174 = cos(pkin(10));
	t181 = r_i_i_C(1) * t174 - r_i_i_C(2) * t172 + pkin(4);
	t188 = r_i_i_C(3) + qJ(5);
	t192 = -(t168 * t188 + t169 * t181) * qJD(4) + t169 * qJD(5);
	t191 = qJD(4) * t181 - qJD(5);
	t189 = t181 * t168 - t188 * t169 + sin(pkin(9)) * pkin(3) + qJ(2);
	t176 = sin(qJ(1));
	t187 = qJD(1) * t176;
	t177 = cos(qJ(1));
	t170 = qJD(1) * t177;
	t186 = qJD(4) * t176;
	t185 = qJD(4) * t177;
	t182 = qJD(1) * t188;
	t180 = -r_i_i_C(1) * t172 - r_i_i_C(2) * t174 - pkin(1) - pkin(7) - qJ(3);
	t179 = qJD(1) * t181;
	t178 = qJD(2) - t192;
	t1 = [-t176 * qJD(3) + t178 * t177 + (-t176 * t189 + t180 * t177) * qJD(1), t170, -t187, (t177 * t179 + t186 * t188) * t169 + (-t176 * t191 + t177 * t182) * t168, t168 * t186 - t169 * t170, 0; t177 * qJD(3) + t178 * t176 + (t180 * t176 + t177 * t189) * qJD(1), t187, t170, (t176 * t179 - t185 * t188) * t169 + (t176 * t182 + t177 * t191) * t168, -t168 * t185 - t169 * t187, 0; 0, 0, 0, t192, qJD(4) * t169, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:44:25
	% EndTime: 2019-10-09 23:44:26
	% DurationCPUTime: 0.27s
	% Computational Cost: add. (262->51), mult. (321->78), div. (0->0), fcn. (256->9), ass. (0->39)
	t218 = pkin(9) + qJ(4);
	t213 = sin(t218);
	t211 = cos(pkin(10)) * pkin(5) + pkin(4);
	t217 = pkin(10) + qJ(6);
	t212 = sin(t217);
	t214 = cos(t217);
	t229 = r_i_i_C(1) * t214 - r_i_i_C(2) * t212 + t211;
	t215 = cos(t218);
	t248 = r_i_i_C(3) + pkin(8) + qJ(5);
	t251 = t248 * t215;
	t257 = -(-t229 * t213 + t251) * qJD(4) - qJD(5) * t213;
	t237 = t248 * t213;
	t256 = t229 * t215 + t237;
	t255 = t251 - pkin(3) * sin(pkin(9)) - t211 * t213 - qJ(2);
	t223 = sin(qJ(1));
	t234 = qJD(1) * t213 + qJD(6);
	t224 = cos(qJ(1));
	t244 = qJD(4) * t224;
	t250 = -t215 * t244 + t234 * t223;
	t245 = qJD(4) * t223;
	t249 = t215 * t245 + t234 * t224;
	t247 = qJD(1) * t223;
	t216 = qJD(1) * t224;
	t246 = qJD(4) * t213;
	t242 = qJD(6) * t215;
	t241 = t215 * qJD(5);
	t235 = -qJD(6) * t213 - qJD(1);
	t233 = -pkin(5) * sin(pkin(10)) - pkin(1) - pkin(7) - qJ(3);
	t232 = r_i_i_C(1) * t212 + r_i_i_C(2) * t214;
	t231 = t235 * t223;
	t230 = t235 * t224;
	t227 = qJD(6) * t232;
	t226 = -t241 + qJD(2) + (t211 * t215 + t237) * qJD(4);
	t225 = qJD(1) * t256;
	t210 = t212 * t231 + t249 * t214;
	t209 = -t249 * t212 + t214 * t231;
	t208 = t212 * t230 - t250 * t214;
	t207 = t250 * t212 + t214 * t230;
	t1 = [t208 * r_i_i_C(1) + t207 * r_i_i_C(2) - t223 * qJD(3) + t226 * t224 + (t255 * t223 + t233 * t224) * qJD(1), t216, -t247, t224 * t225 + (-t232 * t242 - t257) * t223, t213 * t245 - t215 * t216, t209 * r_i_i_C(1) - t210 * r_i_i_C(2); t210 * r_i_i_C(1) + t209 * r_i_i_C(2) + t224 * qJD(3) + t226 * t223 + (t233 * t223 - t255 * t224) * qJD(1), t247, t216, t223 * t225 + (t215 * t227 + t257) * t224, -t213 * t244 - t215 * t247, -t207 * r_i_i_C(1) + t208 * r_i_i_C(2); 0, 0, 0, -t256 * qJD(4) + t213 * t227 + t241, qJD(4) * t215, (t212 * t242 + t214 * t246) * r_i_i_C(2) + (t212 * t246 - t214 * t242) * r_i_i_C(1);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end