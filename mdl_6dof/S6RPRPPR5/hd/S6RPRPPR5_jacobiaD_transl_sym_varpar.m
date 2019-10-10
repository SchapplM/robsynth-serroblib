% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPRPPR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta5]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:22
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRPPR5_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR5_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR5_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPPR5_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPPR5_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR5_jacobiaD_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:22:17
	% EndTime: 2019-10-10 00:22:17
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:22:17
	% EndTime: 2019-10-10 00:22:17
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
	% StartTime: 2019-10-10 00:22:17
	% EndTime: 2019-10-10 00:22:17
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (10->7), mult. (28->12), div. (0->0), fcn. (18->4), ass. (0->5)
	t16 = r_i_i_C(3) + qJ(2);
	t15 = -r_i_i_C(1) * cos(pkin(9)) + r_i_i_C(2) * sin(pkin(9)) - pkin(1);
	t14 = cos(qJ(1));
	t13 = sin(qJ(1));
	t1 = [t14 * qJD(2) + (-t16 * t13 + t15 * t14) * qJD(1), qJD(1) * t14, 0, 0, 0, 0; t13 * qJD(2) + (t15 * t13 + t16 * t14) * qJD(1), qJD(1) * t13, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:22:18
	% EndTime: 2019-10-10 00:22:18
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (43->19), mult. (70->32), div. (0->0), fcn. (46->5), ass. (0->14)
	t35 = r_i_i_C(3) + pkin(7) + qJ(2);
	t26 = sin(qJ(1));
	t34 = qJD(1) * t26;
	t27 = cos(qJ(1));
	t33 = qJD(1) * t27;
	t32 = qJD(3) * t26;
	t31 = qJD(3) * t27;
	t24 = pkin(9) + qJ(3);
	t22 = sin(t24);
	t23 = cos(t24);
	t30 = r_i_i_C(1) * t22 + r_i_i_C(2) * t23;
	t29 = -r_i_i_C(1) * t23 + r_i_i_C(2) * t22 - cos(pkin(9)) * pkin(2) - pkin(1);
	t28 = t30 * qJD(3);
	t1 = [t27 * qJD(2) + t30 * t32 + (-t35 * t26 + t29 * t27) * qJD(1), t33, (t22 * t31 + t23 * t34) * r_i_i_C(2) + (t22 * t34 - t23 * t31) * r_i_i_C(1), 0, 0, 0; t26 * qJD(2) - t27 * t28 + (t29 * t26 + t35 * t27) * qJD(1), t34, (t22 * t32 - t23 * t33) * r_i_i_C(2) + (-t22 * t33 - t23 * t32) * r_i_i_C(1), 0, 0, 0; 0, 0, -t28, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:22:18
	% EndTime: 2019-10-10 00:22:18
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (96->25), mult. (140->39), div. (0->0), fcn. (98->5), ass. (0->16)
	t146 = pkin(9) + qJ(3);
	t144 = sin(t146);
	t145 = cos(t146);
	t158 = r_i_i_C(3) + qJ(4);
	t160 = pkin(3) - r_i_i_C(2);
	t162 = (t160 * t144 - t158 * t145) * qJD(3) - t144 * qJD(4);
	t159 = r_i_i_C(1) + pkin(7) + qJ(2);
	t148 = sin(qJ(1));
	t157 = qJD(1) * t148;
	t149 = cos(qJ(1));
	t156 = qJD(1) * t149;
	t155 = qJD(3) * t145;
	t153 = qJD(3) * t158;
	t152 = -t160 * qJD(3) + qJD(4);
	t151 = -t158 * t144 - t160 * t145 - cos(pkin(9)) * pkin(2) - pkin(1);
	t1 = [t149 * qJD(2) + t162 * t148 + (-t159 * t148 + t151 * t149) * qJD(1), t156, (-t149 * t153 + t160 * t157) * t144 + (t152 * t149 - t158 * t157) * t145, -t144 * t157 + t149 * t155, 0, 0; t148 * qJD(2) - t162 * t149 + (t151 * t148 + t159 * t149) * qJD(1), t157, (-t148 * t153 - t160 * t156) * t144 + (t152 * t148 + t158 * t156) * t145, t144 * t156 + t148 * t155, 0, 0; 0, 0, -t162, qJD(3) * t144, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:22:19
	% EndTime: 2019-10-10 00:22:19
	% DurationCPUTime: 0.21s
	% Computational Cost: add. (153->26), mult. (232->39), div. (0->0), fcn. (176->7), ass. (0->20)
	t171 = sin(pkin(10));
	t172 = cos(pkin(10));
	t170 = pkin(9) + qJ(3);
	t168 = sin(t170);
	t169 = cos(t170);
	t183 = r_i_i_C(1) * t171 + r_i_i_C(2) * t172 + qJ(4);
	t185 = pkin(3) + r_i_i_C(3) + qJ(5);
	t180 = t168 * t185 - t169 * t183;
	t176 = -qJD(3) * t180 + t168 * qJD(4) + t169 * qJD(5);
	t191 = t176 + (r_i_i_C(1) * t172 - r_i_i_C(2) * t171 + pkin(4) + pkin(7) + qJ(2)) * qJD(1);
	t174 = sin(qJ(1));
	t189 = qJD(1) * t174;
	t175 = cos(qJ(1));
	t188 = qJD(1) * t175;
	t187 = qJD(3) * t174;
	t186 = qJD(3) * t175;
	t181 = -t168 * t183 - t169 * t185;
	t178 = qJD(2) + (-cos(pkin(9)) * pkin(2) - pkin(1) + t181) * qJD(1);
	t177 = qJD(3) * t181 + qJD(4) * t169 - qJD(5) * t168;
	t1 = [-t174 * t191 + t178 * t175, t188, t175 * t177 + t180 * t189, -t168 * t189 + t169 * t186, -t168 * t186 - t169 * t189, 0; t178 * t174 + t175 * t191, t189, t174 * t177 - t180 * t188, t168 * t188 + t169 * t187, -t168 * t187 + t169 * t188, 0; 0, 0, t176, qJD(3) * t168, qJD(3) * t169, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:22:19
	% EndTime: 2019-10-10 00:22:19
	% DurationCPUTime: 0.27s
	% Computational Cost: add. (303->50), mult. (379->75), div. (0->0), fcn. (300->9), ass. (0->38)
	t222 = pkin(9) + qJ(3);
	t218 = sin(t222);
	t220 = cos(t222);
	t241 = pkin(5) * sin(pkin(10)) + qJ(4);
	t245 = t220 * qJD(5);
	t244 = pkin(3) + r_i_i_C(3) + pkin(8) + qJ(5);
	t254 = t244 * t218;
	t261 = (-t241 * t220 + t254) * qJD(3) - (cos(pkin(10)) * pkin(5) + pkin(4) + pkin(7) + qJ(2)) * qJD(1) - t218 * qJD(4) - t245;
	t221 = pkin(10) + qJ(6);
	t217 = sin(t221);
	t219 = cos(t221);
	t234 = r_i_i_C(1) * t217 + r_i_i_C(2) * t219 + t241;
	t259 = -t234 * t220 + t254;
	t226 = sin(qJ(1));
	t240 = qJD(6) * t218 + qJD(1);
	t258 = t226 * t240;
	t227 = cos(qJ(1));
	t257 = t227 * t240;
	t251 = qJD(1) * t226;
	t250 = qJD(1) * t227;
	t249 = qJD(3) * t218;
	t248 = qJD(3) * t226;
	t247 = qJD(3) * t227;
	t246 = qJD(6) * t220;
	t243 = t220 * t248;
	t242 = t220 * t247;
	t239 = -qJD(1) * t218 - qJD(6);
	t237 = t244 * t220;
	t233 = qJD(4) + (r_i_i_C(1) * t219 - r_i_i_C(2) * t217) * qJD(6);
	t232 = t239 * t227 - t243;
	t231 = t239 * t226 + t242;
	t230 = qJD(2) + (-t241 * t218 - cos(pkin(9)) * pkin(2) - pkin(1) - t237) * qJD(1);
	t228 = -qJD(5) * t218 + t233 * t220 + (-t234 * t218 - t237) * qJD(3);
	t214 = t231 * t217 + t219 * t257;
	t213 = -t217 * t257 + t231 * t219;
	t212 = t232 * t217 - t219 * t258;
	t211 = t217 * t258 + t232 * t219;
	t1 = [t212 * r_i_i_C(1) + t211 * r_i_i_C(2) + t261 * t226 + t230 * t227, t250, t228 * t227 + t259 * t251, -t218 * t251 + t242, -t218 * t247 - t220 * t251, t213 * r_i_i_C(1) - t214 * r_i_i_C(2); t214 * r_i_i_C(1) + t213 * r_i_i_C(2) + t230 * t226 - t261 * t227, t251, t228 * t226 - t250 * t259, t218 * t250 + t243, -t218 * t248 + t220 * t250, -t211 * r_i_i_C(1) + t212 * r_i_i_C(2); 0, 0, -qJD(3) * t259 + t233 * t218 + t245, t249, qJD(3) * t220, (-t217 * t249 + t219 * t246) * r_i_i_C(2) + (t217 * t246 + t219 * t249) * r_i_i_C(1);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end