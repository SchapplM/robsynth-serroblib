% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RPRPR13
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RPRPR13_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR13_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR13_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RPRPR13_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRPR13_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR13_jacobiaD_transl_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:33:40
	% EndTime: 2019-12-31 18:33:40
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:33:40
	% EndTime: 2019-12-31 18:33:40
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:33:40
	% EndTime: 2019-12-31 18:33:40
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (10->7), mult. (28->12), div. (0->0), fcn. (18->4), ass. (0->5)
	t16 = r_i_i_C(3) + qJ(2);
	t15 = -r_i_i_C(1) * cos(pkin(8)) + r_i_i_C(2) * sin(pkin(8)) - pkin(1);
	t14 = cos(qJ(1));
	t13 = sin(qJ(1));
	t1 = [t14 * qJD(2) + (-t16 * t13 + t15 * t14) * qJD(1), qJD(1) * t14, 0, 0, 0; t13 * qJD(2) + (t15 * t13 + t16 * t14) * qJD(1), qJD(1) * t13, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:33:40
	% EndTime: 2019-12-31 18:33:40
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (43->19), mult. (70->32), div. (0->0), fcn. (46->5), ass. (0->14)
	t35 = r_i_i_C(3) + pkin(6) + qJ(2);
	t26 = sin(qJ(1));
	t34 = qJD(1) * t26;
	t27 = cos(qJ(1));
	t33 = qJD(1) * t27;
	t32 = qJD(3) * t26;
	t31 = qJD(3) * t27;
	t24 = pkin(8) + qJ(3);
	t22 = sin(t24);
	t23 = cos(t24);
	t30 = r_i_i_C(1) * t22 + r_i_i_C(2) * t23;
	t29 = -r_i_i_C(1) * t23 + r_i_i_C(2) * t22 - cos(pkin(8)) * pkin(2) - pkin(1);
	t28 = t30 * qJD(3);
	t1 = [t27 * qJD(2) + t30 * t32 + (-t35 * t26 + t29 * t27) * qJD(1), t33, (t22 * t31 + t23 * t34) * r_i_i_C(2) + (t22 * t34 - t23 * t31) * r_i_i_C(1), 0, 0; t26 * qJD(2) - t27 * t28 + (t29 * t26 + t35 * t27) * qJD(1), t34, (t22 * t32 - t23 * t33) * r_i_i_C(2) + (-t22 * t33 - t23 * t32) * r_i_i_C(1), 0, 0; 0, 0, -t28, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:33:41
	% EndTime: 2019-12-31 18:33:41
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (96->25), mult. (140->39), div. (0->0), fcn. (98->5), ass. (0->16)
	t146 = pkin(8) + qJ(3);
	t144 = sin(t146);
	t145 = cos(t146);
	t158 = r_i_i_C(3) + qJ(4);
	t160 = pkin(3) - r_i_i_C(2);
	t162 = (t160 * t144 - t158 * t145) * qJD(3) - t144 * qJD(4);
	t159 = r_i_i_C(1) + pkin(6) + qJ(2);
	t148 = sin(qJ(1));
	t157 = qJD(1) * t148;
	t149 = cos(qJ(1));
	t156 = qJD(1) * t149;
	t155 = qJD(3) * t145;
	t153 = qJD(3) * t158;
	t152 = -t160 * qJD(3) + qJD(4);
	t151 = -t158 * t144 - t160 * t145 - cos(pkin(8)) * pkin(2) - pkin(1);
	t1 = [t149 * qJD(2) + t162 * t148 + (-t159 * t148 + t151 * t149) * qJD(1), t156, (-t149 * t153 + t160 * t157) * t144 + (t152 * t149 - t158 * t157) * t145, -t144 * t157 + t149 * t155, 0; t148 * qJD(2) - t162 * t149 + (t151 * t148 + t159 * t149) * qJD(1), t157, (-t148 * t153 - t160 * t156) * t144 + (t152 * t148 + t158 * t156) * t145, t144 * t156 + t148 * t155, 0; 0, 0, -t162, qJD(3) * t144, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:33:41
	% EndTime: 2019-12-31 18:33:41
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (197->47), mult. (324->78), div. (0->0), fcn. (254->7), ass. (0->35)
	t205 = pkin(8) + qJ(3);
	t203 = sin(t205);
	t204 = cos(t205);
	t223 = pkin(3) + pkin(7) + r_i_i_C(3);
	t220 = t223 * t203;
	t236 = (-qJ(4) * t204 + t220) * qJD(3) - t203 * qJD(4);
	t209 = cos(qJ(5));
	t219 = qJD(5) * t203 + qJD(1);
	t234 = t209 * t219;
	t207 = sin(qJ(5));
	t233 = t219 * t207;
	t232 = pkin(4) + pkin(6) + qJ(2);
	t208 = sin(qJ(1));
	t230 = qJD(1) * t208;
	t210 = cos(qJ(1));
	t229 = qJD(1) * t210;
	t228 = qJD(3) * t203;
	t227 = qJD(3) * t204;
	t226 = qJD(3) * t209;
	t225 = qJD(5) * t204;
	t222 = t208 * t227;
	t221 = t210 * t227;
	t218 = -qJD(1) * t203 - qJD(5);
	t217 = t218 * t210;
	t216 = r_i_i_C(1) * t207 + r_i_i_C(2) * t209 + qJ(4);
	t215 = qJD(3) * t216;
	t214 = -qJ(4) * t203 - t223 * t204 - cos(pkin(8)) * pkin(2) - pkin(1);
	t213 = qJD(4) + (r_i_i_C(1) * t209 - r_i_i_C(2) * t207) * qJD(5);
	t212 = t218 * t208 + t221;
	t211 = -t223 * qJD(3) + t213;
	t201 = t212 * t207 + t210 * t234;
	t200 = t212 * t209 - t210 * t233;
	t199 = -t208 * t234 + (t217 - t222) * t207;
	t198 = t209 * t217 + (-t204 * t226 + t233) * t208;
	t1 = [t199 * r_i_i_C(1) + t198 * r_i_i_C(2) + t210 * qJD(2) + t236 * t208 + (-t232 * t208 + t214 * t210) * qJD(1), t229, (-t210 * t215 + t223 * t230) * t203 + (t211 * t210 - t216 * t230) * t204, -t203 * t230 + t221, t200 * r_i_i_C(1) - t201 * r_i_i_C(2); t201 * r_i_i_C(1) + t200 * r_i_i_C(2) + t208 * qJD(2) - t236 * t210 + (t214 * t208 + t232 * t210) * qJD(1), t230, (-t208 * t215 - t223 * t229) * t203 + (t211 * t208 + t216 * t229) * t204, t203 * t229 + t222, -t198 * r_i_i_C(1) + t199 * r_i_i_C(2); 0, 0, t213 * t203 + (t216 * t204 - t220) * qJD(3), t228, (-t207 * t228 + t209 * t225) * r_i_i_C(2) + (t203 * t226 + t207 * t225) * r_i_i_C(1);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end