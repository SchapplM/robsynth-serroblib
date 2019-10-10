% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPRPPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4,theta5]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:15
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRPPR1_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR1_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR1_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPPR1_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPPR1_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPPR1_jacobiaD_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:15:13
	% EndTime: 2019-10-10 00:15:13
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:15:13
	% EndTime: 2019-10-10 00:15:13
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
	% StartTime: 2019-10-10 00:15:13
	% EndTime: 2019-10-10 00:15:13
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
	% StartTime: 2019-10-10 00:15:13
	% EndTime: 2019-10-10 00:15:13
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (43->18), mult. (68->31), div. (0->0), fcn. (42->6), ass. (0->13)
	t24 = sin(qJ(3));
	t25 = cos(qJ(3));
	t26 = (r_i_i_C(1) * t24 + r_i_i_C(2) * t25) * qJD(3);
	t33 = pkin(7) + r_i_i_C(3);
	t32 = qJD(1) * t24;
	t31 = qJD(1) * t25;
	t30 = qJD(3) * t24;
	t29 = qJD(3) * t25;
	t27 = -r_i_i_C(1) * t25 + r_i_i_C(2) * t24 - pkin(2);
	t23 = qJ(1) + pkin(9);
	t22 = cos(t23);
	t21 = sin(t23);
	t1 = [t21 * t26 + (-cos(qJ(1)) * pkin(1) - t33 * t21 + t27 * t22) * qJD(1), 0, (t21 * t31 + t22 * t30) * r_i_i_C(2) + (t21 * t32 - t22 * t29) * r_i_i_C(1), 0, 0, 0; -t22 * t26 + (-sin(qJ(1)) * pkin(1) + t33 * t22 + t27 * t21) * qJD(1), 0, (t21 * t30 - t22 * t31) * r_i_i_C(2) + (-t21 * t29 - t22 * t32) * r_i_i_C(1), 0, 0, 0; 0, 0, -t26, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:15:14
	% EndTime: 2019-10-10 00:15:14
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (84->19), mult. (94->27), div. (0->0), fcn. (61->8), ass. (0->15)
	t33 = qJ(3) + pkin(10);
	t29 = sin(t33);
	t31 = cos(t33);
	t47 = -r_i_i_C(1) * t31 + r_i_i_C(2) * t29 - cos(qJ(3)) * pkin(3);
	t45 = r_i_i_C(3) + qJ(4) + pkin(7);
	t34 = qJ(1) + pkin(9);
	t32 = cos(t34);
	t44 = qJD(1) * t32;
	t42 = -pkin(2) + t47;
	t41 = sin(qJ(3)) * pkin(3) + r_i_i_C(1) * t29 + r_i_i_C(2) * t31;
	t30 = sin(t34);
	t40 = t41 * t30;
	t39 = qJD(3) * t47;
	t38 = t41 * qJD(3);
	t1 = [t32 * qJD(4) + qJD(3) * t40 + (-cos(qJ(1)) * pkin(1) - t45 * t30 + t42 * t32) * qJD(1), 0, qJD(1) * t40 + t32 * t39, t44, 0, 0; t30 * qJD(4) - t32 * t38 + (-sin(qJ(1)) * pkin(1) + t45 * t32 + t42 * t30) * qJD(1), 0, t30 * t39 - t41 * t44, qJD(1) * t30, 0, 0; 0, 0, -t38, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:15:15
	% EndTime: 2019-10-10 00:15:15
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (191->28), mult. (208->38), div. (0->0), fcn. (154->10), ass. (0->20)
	t186 = qJ(3) + pkin(10);
	t182 = sin(t186);
	t184 = cos(t186);
	t188 = sin(pkin(11));
	t189 = cos(pkin(11));
	t200 = r_i_i_C(1) * t189 - r_i_i_C(2) * t188 + pkin(4);
	t205 = r_i_i_C(3) + qJ(5);
	t195 = t200 * t182 - t205 * t184 + sin(qJ(3)) * pkin(3);
	t209 = t195 * qJD(3) - t182 * qJD(5);
	t208 = -t205 * t182 - t200 * t184 - cos(qJ(3)) * pkin(3);
	t187 = qJ(1) + pkin(9);
	t183 = sin(t187);
	t204 = qJD(1) * t183;
	t185 = cos(t187);
	t203 = qJD(1) * t185;
	t202 = qJD(3) * t184;
	t199 = t188 * r_i_i_C(1) + t189 * r_i_i_C(2) + pkin(7) + qJ(4);
	t197 = -pkin(2) + t208;
	t194 = t208 * qJD(3) + qJD(5) * t184;
	t1 = [t185 * qJD(4) + t209 * t183 + (-cos(qJ(1)) * pkin(1) - t199 * t183 + t197 * t185) * qJD(1), 0, t194 * t185 + t195 * t204, t203, -t182 * t204 + t185 * t202, 0; t183 * qJD(4) - t209 * t185 + (-sin(qJ(1)) * pkin(1) + t199 * t185 + t197 * t183) * qJD(1), 0, t194 * t183 - t195 * t203, t204, t182 * t203 + t183 * t202, 0; 0, 0, -t209, 0, qJD(3) * t182, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:15:14
	% EndTime: 2019-10-10 00:15:15
	% DurationCPUTime: 0.29s
	% Computational Cost: add. (366->52), mult. (333->79), div. (0->0), fcn. (263->12), ass. (0->38)
	t221 = cos(pkin(11)) * pkin(5) + pkin(4);
	t230 = qJ(3) + pkin(10);
	t224 = sin(t230);
	t227 = cos(t230);
	t260 = r_i_i_C(3) + pkin(8) + qJ(5);
	t242 = t260 * t227 - sin(qJ(3)) * pkin(3);
	t250 = t224 * qJD(5);
	t266 = t250 + (-t221 * t224 + t242) * qJD(3);
	t229 = pkin(11) + qJ(6);
	t223 = sin(t229);
	t226 = cos(t229);
	t243 = r_i_i_C(1) * t226 - r_i_i_C(2) * t223 + t221;
	t238 = -t243 * t224 + t242;
	t263 = -t260 * t224 - cos(qJ(3)) * pkin(3);
	t231 = qJ(1) + pkin(9);
	t228 = cos(t231);
	t258 = t226 * t228;
	t225 = sin(t231);
	t257 = qJD(1) * t225;
	t256 = qJD(1) * t228;
	t255 = qJD(3) * t224;
	t254 = qJD(3) * t227;
	t253 = qJD(3) * t228;
	t252 = qJD(6) * t224;
	t251 = qJD(6) * t227;
	t248 = pkin(5) * sin(pkin(11)) + qJ(4) + pkin(7);
	t247 = -qJD(1) + t251;
	t246 = qJD(1) * t227 - qJD(6);
	t245 = r_i_i_C(1) * t223 + r_i_i_C(2) * t226;
	t244 = t247 * t223;
	t240 = -t221 * t227 - pkin(2) + t263;
	t239 = t224 * t253 + t246 * t225;
	t237 = qJD(5) * t227 + t245 * t252 + (-t243 * t227 + t263) * qJD(3);
	t220 = -t246 * t258 + (t226 * t255 + t244) * t225;
	t219 = t247 * t226 * t225 + (-t225 * t255 + t246 * t228) * t223;
	t218 = t239 * t226 + t228 * t244;
	t217 = t239 * t223 - t247 * t258;
	t1 = [t220 * r_i_i_C(1) + t219 * r_i_i_C(2) + t228 * qJD(4) - t266 * t225 + (-cos(qJ(1)) * pkin(1) - t248 * t225 + t240 * t228) * qJD(1), 0, t237 * t228 - t238 * t257, t256, -t224 * t257 + t227 * t253, t217 * r_i_i_C(1) + t218 * r_i_i_C(2); -t218 * r_i_i_C(1) + t217 * r_i_i_C(2) + t225 * qJD(4) + t266 * t228 + (-sin(qJ(1)) * pkin(1) + t248 * t228 + t240 * t225) * qJD(1), 0, t237 * t225 + t238 * t256, t257, t224 * t256 + t225 * t254, -t219 * r_i_i_C(1) + t220 * r_i_i_C(2); 0, 0, t238 * qJD(3) - t245 * t251 + t250, 0, t255, (t223 * t252 - t226 * t254) * r_i_i_C(2) + (-t223 * t254 - t226 * t252) * r_i_i_C(1);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end