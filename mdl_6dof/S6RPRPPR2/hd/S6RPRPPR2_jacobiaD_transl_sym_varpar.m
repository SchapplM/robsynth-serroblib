% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPRPPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:17
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRPPR2_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR2_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR2_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPPR2_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPPR2_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR2_jacobiaD_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:17:00
	% EndTime: 2019-10-10 00:17:00
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:17:00
	% EndTime: 2019-10-10 00:17:00
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
	% StartTime: 2019-10-10 00:17:00
	% EndTime: 2019-10-10 00:17:00
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
	% StartTime: 2019-10-10 00:17:00
	% EndTime: 2019-10-10 00:17:00
	% DurationCPUTime: 0.10s
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
	% StartTime: 2019-10-10 00:17:00
	% EndTime: 2019-10-10 00:17:00
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
	% StartTime: 2019-10-10 00:17:01
	% EndTime: 2019-10-10 00:17:01
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (161->25), mult. (164->34), div. (0->0), fcn. (113->8), ass. (0->18)
	t159 = qJ(3) + pkin(10);
	t155 = sin(t159);
	t157 = cos(t159);
	t174 = r_i_i_C(3) + qJ(5);
	t178 = pkin(4) - r_i_i_C(2);
	t166 = t178 * t155 - t174 * t157 + sin(qJ(3)) * pkin(3);
	t180 = t166 * qJD(3) - t155 * qJD(5);
	t179 = -t174 * t155 - t178 * t157 - cos(qJ(3)) * pkin(3);
	t175 = r_i_i_C(1) + qJ(4) + pkin(7);
	t160 = qJ(1) + pkin(9);
	t156 = sin(t160);
	t173 = qJD(1) * t156;
	t158 = cos(t160);
	t172 = qJD(1) * t158;
	t171 = qJD(3) * t157;
	t168 = -pkin(2) + t179;
	t165 = t179 * qJD(3) + qJD(5) * t157;
	t1 = [t158 * qJD(4) + t180 * t156 + (-cos(qJ(1)) * pkin(1) - t175 * t156 + t168 * t158) * qJD(1), 0, t165 * t158 + t166 * t173, t172, -t155 * t173 + t158 * t171, 0; t156 * qJD(4) - t180 * t158 + (-sin(qJ(1)) * pkin(1) + t175 * t158 + t168 * t156) * qJD(1), 0, t165 * t156 - t166 * t172, t173, t155 * t172 + t156 * t171, 0; 0, 0, -t180, 0, qJD(3) * t155, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:17:01
	% EndTime: 2019-10-10 00:17:01
	% DurationCPUTime: 0.28s
	% Computational Cost: add. (314->47), mult. (348->76), div. (0->0), fcn. (269->10), ass. (0->35)
	t219 = qJ(3) + pkin(10);
	t215 = sin(t219);
	t217 = cos(t219);
	t240 = pkin(4) + pkin(8) + r_i_i_C(3);
	t231 = t240 * t215 + sin(qJ(3)) * pkin(3);
	t258 = -t215 * qJD(5) + (-qJ(5) * t217 + t231) * qJD(3);
	t222 = sin(qJ(6));
	t224 = cos(qJ(6));
	t232 = r_i_i_C(1) * t222 + r_i_i_C(2) * t224 + qJ(5);
	t256 = -t232 * t217 + t231;
	t254 = -t240 * t217 - cos(qJ(3)) * pkin(3);
	t236 = qJD(6) * t215 + qJD(1);
	t243 = qJD(3) * t224;
	t253 = -t217 * t243 + t236 * t222;
	t244 = qJD(3) * t222;
	t252 = t217 * t244 + t236 * t224;
	t249 = pkin(5) + qJ(4) + pkin(7);
	t220 = qJ(1) + pkin(9);
	t216 = sin(t220);
	t247 = qJD(1) * t216;
	t218 = cos(t220);
	t246 = qJD(1) * t218;
	t245 = qJD(3) * t217;
	t242 = qJD(6) * t217;
	t235 = -qJD(1) * t215 - qJD(6);
	t234 = t235 * t222;
	t233 = t235 * t224;
	t229 = -qJ(5) * t215 - pkin(2) + t254;
	t228 = qJD(5) + (r_i_i_C(1) * t224 - r_i_i_C(2) * t222) * qJD(6);
	t226 = t228 * t217 + (-t232 * t215 + t254) * qJD(3);
	t213 = t216 * t234 + t252 * t218;
	t212 = t216 * t233 - t253 * t218;
	t211 = -t252 * t216 + t218 * t234;
	t210 = t253 * t216 + t218 * t233;
	t1 = [t211 * r_i_i_C(1) + t210 * r_i_i_C(2) + t218 * qJD(4) + t258 * t216 + (-cos(qJ(1)) * pkin(1) - t249 * t216 + t229 * t218) * qJD(1), 0, t226 * t218 + t256 * t247, t246, -t215 * t247 + t218 * t245, t212 * r_i_i_C(1) - t213 * r_i_i_C(2); t213 * r_i_i_C(1) + t212 * r_i_i_C(2) + t216 * qJD(4) - t258 * t218 + (-sin(qJ(1)) * pkin(1) + t249 * t218 + t229 * t216) * qJD(1), 0, t226 * t216 - t246 * t256, t247, t215 * t246 + t216 * t245, -t210 * r_i_i_C(1) + t211 * r_i_i_C(2); 0, 0, -qJD(3) * t256 + t228 * t215, 0, qJD(3) * t215, (-t215 * t244 + t224 * t242) * r_i_i_C(2) + (t215 * t243 + t222 * t242) * r_i_i_C(1);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end