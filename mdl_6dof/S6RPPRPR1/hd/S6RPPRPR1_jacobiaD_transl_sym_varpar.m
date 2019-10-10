% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPPRPR1
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3,theta5]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:34
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPPRPR1_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR1_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR1_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRPR1_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRPR1_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRPR1_jacobiaD_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:34:17
	% EndTime: 2019-10-09 23:34:17
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:34:17
	% EndTime: 2019-10-09 23:34:17
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
	% StartTime: 2019-10-09 23:34:17
	% EndTime: 2019-10-09 23:34:17
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
	% StartTime: 2019-10-09 23:34:17
	% EndTime: 2019-10-09 23:34:17
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (26->10), mult. (32->14), div. (0->0), fcn. (20->6), ass. (0->6)
	t19 = r_i_i_C(3) + qJ(3);
	t18 = -r_i_i_C(1) * cos(pkin(10)) + r_i_i_C(2) * sin(pkin(10)) - pkin(2);
	t15 = qJ(1) + pkin(9);
	t14 = cos(t15);
	t13 = sin(t15);
	t1 = [t14 * qJD(3) + (-cos(qJ(1)) * pkin(1) - t19 * t13 + t18 * t14) * qJD(1), 0, qJD(1) * t14, 0, 0, 0; t13 * qJD(3) + (-sin(qJ(1)) * pkin(1) + t19 * t14 + t18 * t13) * qJD(1), 0, qJD(1) * t13, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:34:17
	% EndTime: 2019-10-09 23:34:17
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (71->22), mult. (74->34), div. (0->0), fcn. (48->7), ass. (0->15)
	t40 = r_i_i_C(3) + pkin(7) + qJ(3);
	t31 = qJ(1) + pkin(9);
	t27 = sin(t31);
	t39 = qJD(1) * t27;
	t29 = cos(t31);
	t38 = qJD(1) * t29;
	t37 = qJD(4) * t27;
	t36 = qJD(4) * t29;
	t30 = pkin(10) + qJ(4);
	t26 = sin(t30);
	t28 = cos(t30);
	t35 = r_i_i_C(1) * t26 + r_i_i_C(2) * t28;
	t34 = -r_i_i_C(1) * t28 + r_i_i_C(2) * t26 - cos(pkin(10)) * pkin(3) - pkin(2);
	t33 = t35 * qJD(4);
	t1 = [t29 * qJD(3) + t35 * t37 + (-cos(qJ(1)) * pkin(1) - t40 * t27 + t34 * t29) * qJD(1), 0, t38, (t26 * t36 + t28 * t39) * r_i_i_C(2) + (t26 * t39 - t28 * t36) * r_i_i_C(1), 0, 0; t27 * qJD(3) - t29 * t33 + (-sin(qJ(1)) * pkin(1) + t40 * t29 + t34 * t27) * qJD(1), 0, t39, (t26 * t37 - t28 * t38) * r_i_i_C(2) + (-t26 * t38 - t28 * t37) * r_i_i_C(1), 0, 0; 0, 0, 0, -t33, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:34:19
	% EndTime: 2019-10-09 23:34:19
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (178->30), mult. (188->43), div. (0->0), fcn. (141->9), ass. (0->19)
	t183 = pkin(10) + qJ(4);
	t179 = sin(t183);
	t181 = cos(t183);
	t185 = sin(pkin(11));
	t186 = cos(pkin(11));
	t193 = r_i_i_C(1) * t186 - r_i_i_C(2) * t185 + pkin(4);
	t198 = r_i_i_C(3) + qJ(5);
	t199 = t193 * t179 - t198 * t181;
	t200 = t199 * qJD(4) - t179 * qJD(5);
	t184 = qJ(1) + pkin(9);
	t180 = sin(t184);
	t197 = qJD(1) * t180;
	t182 = cos(t184);
	t196 = qJD(1) * t182;
	t195 = qJD(4) * t182;
	t192 = t185 * r_i_i_C(1) + t186 * r_i_i_C(2) + pkin(7) + qJ(3);
	t191 = -t198 * t179 - t193 * t181;
	t189 = -cos(pkin(10)) * pkin(3) - pkin(2) + t191;
	t1 = [t182 * qJD(3) + t200 * t180 + (-cos(qJ(1)) * pkin(1) - t192 * t180 + t189 * t182) * qJD(1), 0, t196, (t193 * t197 - t198 * t195) * t179 + (-t198 * t197 + (-t193 * qJD(4) + qJD(5)) * t182) * t181, -t179 * t197 + t181 * t195, 0; t180 * qJD(3) - t200 * t182 + (-sin(qJ(1)) * pkin(1) + t192 * t182 + t189 * t180) * qJD(1), 0, t197, -t199 * t196 + (t191 * qJD(4) + qJD(5) * t181) * t180, t180 * qJD(4) * t181 + t179 * t196, 0; 0, 0, 0, -t200, qJD(4) * t179, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:34:19
	% EndTime: 2019-10-09 23:34:19
	% DurationCPUTime: 0.26s
	% Computational Cost: add. (353->51), mult. (313->78), div. (0->0), fcn. (250->11), ass. (0->38)
	t218 = cos(pkin(11)) * pkin(5) + pkin(4);
	t227 = pkin(10) + qJ(4);
	t221 = sin(t227);
	t245 = t221 * qJD(5);
	t224 = cos(t227);
	t255 = r_i_i_C(3) + pkin(8) + qJ(5);
	t256 = t255 * t224;
	t259 = (-t218 * t221 + t256) * qJD(4) + t245;
	t226 = pkin(11) + qJ(6);
	t220 = sin(t226);
	t223 = cos(t226);
	t236 = r_i_i_C(1) * t223 - r_i_i_C(2) * t220 + t218;
	t233 = -t236 * t221 + t256;
	t228 = qJ(1) + pkin(9);
	t225 = cos(t228);
	t253 = t223 * t225;
	t222 = sin(t228);
	t252 = qJD(1) * t222;
	t251 = qJD(1) * t225;
	t250 = qJD(4) * t221;
	t249 = qJD(4) * t224;
	t248 = qJD(4) * t225;
	t247 = qJD(6) * t221;
	t246 = qJD(6) * t224;
	t244 = t255 * t221;
	t241 = pkin(5) * sin(pkin(11)) + pkin(7) + qJ(3);
	t240 = -qJD(1) + t246;
	t239 = qJD(1) * t224 - qJD(6);
	t238 = r_i_i_C(1) * t220 + r_i_i_C(2) * t223;
	t237 = t240 * t220;
	t235 = -t218 * t224 - cos(pkin(10)) * pkin(3) - pkin(2) - t244;
	t234 = t221 * t248 + t239 * t222;
	t232 = qJD(5) * t224 + t238 * t247 + (-t236 * t224 - t244) * qJD(4);
	t217 = -t239 * t253 + (t223 * t250 + t237) * t222;
	t216 = t240 * t223 * t222 + (-t222 * t250 + t239 * t225) * t220;
	t215 = t234 * t223 + t225 * t237;
	t214 = t234 * t220 - t240 * t253;
	t1 = [t217 * r_i_i_C(1) + t216 * r_i_i_C(2) + t225 * qJD(3) - t259 * t222 + (-cos(qJ(1)) * pkin(1) - t241 * t222 + t235 * t225) * qJD(1), 0, t251, t232 * t225 - t233 * t252, -t221 * t252 + t224 * t248, t214 * r_i_i_C(1) + t215 * r_i_i_C(2); -t215 * r_i_i_C(1) + t214 * r_i_i_C(2) + t222 * qJD(3) + t259 * t225 + (-sin(qJ(1)) * pkin(1) + t241 * t225 + t235 * t222) * qJD(1), 0, t252, t232 * t222 + t233 * t251, t221 * t251 + t222 * t249, -t216 * r_i_i_C(1) + t217 * r_i_i_C(2); 0, 0, 0, t233 * qJD(4) - t238 * t246 + t245, t250, (t220 * t247 - t223 * t249) * r_i_i_C(2) + (-t220 * t249 - t223 * t247) * r_i_i_C(1);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end