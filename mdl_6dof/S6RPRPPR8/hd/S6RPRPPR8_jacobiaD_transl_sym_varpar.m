% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPRPPR8
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:27
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRPPR8_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR8_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR8_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPPR8_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPPR8_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPPR8_jacobiaD_transl_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:27:23
	% EndTime: 2019-10-10 00:27:23
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:27:23
	% EndTime: 2019-10-10 00:27:23
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
	% StartTime: 2019-10-10 00:27:23
	% EndTime: 2019-10-10 00:27:23
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
	% StartTime: 2019-10-10 00:27:23
	% EndTime: 2019-10-10 00:27:23
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
	% StartTime: 2019-10-10 00:27:24
	% EndTime: 2019-10-10 00:27:24
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (48->21), mult. (142->35), div. (0->0), fcn. (98->4), ass. (0->17)
	t134 = sin(qJ(3));
	t136 = cos(qJ(3));
	t146 = r_i_i_C(3) + qJ(4);
	t147 = pkin(3) + r_i_i_C(1);
	t151 = -(t146 * t134 + t147 * t136) * qJD(3) + t136 * qJD(4);
	t150 = t147 * qJD(3) - qJD(4);
	t148 = t147 * t134 - t146 * t136 + qJ(2);
	t135 = sin(qJ(1));
	t145 = qJD(1) * t135;
	t137 = cos(qJ(1));
	t144 = qJD(1) * t137;
	t143 = qJD(3) * t135;
	t142 = qJD(3) * t137;
	t140 = qJD(1) * t147;
	t139 = qJD(1) * t146;
	t138 = qJD(2) + (-pkin(1) - pkin(7) - r_i_i_C(2)) * qJD(1) - t151;
	t1 = [t138 * t137 - t148 * t145, t144, (t137 * t140 + t146 * t143) * t136 + (-t150 * t135 + t137 * t139) * t134, t134 * t143 - t136 * t144, 0, 0; t138 * t135 + t148 * t144, t145, (t135 * t140 - t146 * t142) * t136 + (t135 * t139 + t150 * t137) * t134, -t134 * t142 - t136 * t145, 0, 0; 0, 0, t151, qJD(3) * t136, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:27:23
	% EndTime: 2019-10-10 00:27:23
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (62->24), mult. (176->36), div. (0->0), fcn. (121->4), ass. (0->17)
	t20 = sin(qJ(3));
	t22 = cos(qJ(3));
	t27 = pkin(3) + pkin(4) - r_i_i_C(2);
	t32 = r_i_i_C(1) + qJ(4);
	t37 = -(t32 * t20 + t27 * t22) * qJD(3) + t22 * qJD(4);
	t36 = -qJD(5) + (t27 * t20 - t32 * t22 + qJ(2)) * qJD(1);
	t35 = t27 * qJD(3) - qJD(4);
	t21 = sin(qJ(1));
	t19 = qJD(1) * t21;
	t23 = cos(qJ(1));
	t31 = qJD(1) * t23;
	t30 = qJD(3) * t21;
	t29 = qJD(3) * t23;
	t26 = qJD(1) * t32;
	t25 = qJD(1) * t27;
	t24 = qJD(2) + (-pkin(1) - pkin(7) + r_i_i_C(3) + qJ(5)) * qJD(1) - t37;
	t1 = [-t36 * t21 + t24 * t23, t31, (t23 * t25 + t32 * t30) * t22 + (-t35 * t21 + t23 * t26) * t20, t20 * t30 - t22 * t31, t19, 0; t24 * t21 + t36 * t23, t19, (t21 * t25 - t32 * t29) * t22 + (t21 * t26 + t35 * t23) * t20, -t22 * t19 - t20 * t29, -t31, 0; 0, 0, t37, qJD(3) * t22, 0, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:27:24
	% EndTime: 2019-10-10 00:27:24
	% DurationCPUTime: 0.30s
	% Computational Cost: add. (126->47), mult. (382->75), div. (0->0), fcn. (292->6), ass. (0->36)
	t201 = sin(qJ(3));
	t204 = cos(qJ(3));
	t218 = pkin(3) + pkin(4) + pkin(8) + r_i_i_C(3);
	t225 = pkin(5) + qJ(4);
	t231 = -qJD(5) + (t218 * t201 - t225 * t204 + qJ(2)) * qJD(1);
	t203 = cos(qJ(6));
	t215 = qJD(6) * t204 + qJD(1);
	t230 = t203 * t215;
	t200 = sin(qJ(6));
	t229 = t215 * t200;
	t228 = t218 * t204;
	t208 = -(r_i_i_C(1) * t200 + r_i_i_C(2) * t203) * qJD(6) + qJD(4);
	t227 = -t218 * qJD(3) + t208;
	t224 = qJD(1) * t204;
	t205 = cos(qJ(1));
	t223 = qJD(1) * t205;
	t202 = sin(qJ(1));
	t222 = qJD(3) * t202;
	t221 = qJD(3) * t204;
	t220 = qJD(3) * t205;
	t219 = qJD(6) * t201;
	t217 = t201 * t222;
	t216 = t201 * t220;
	t214 = qJD(6) + t224;
	t213 = qJD(1) * t218;
	t211 = t214 * t205;
	t210 = r_i_i_C(1) * t203 - r_i_i_C(2) * t200 + t225;
	t209 = qJD(1) * t210;
	t207 = t214 * t202 + t216;
	t206 = -t204 * qJD(4) + qJD(2) + (-pkin(1) - pkin(7) + qJ(5)) * qJD(1) + (t225 * t201 + t228) * qJD(3);
	t199 = qJD(1) * t202;
	t198 = t203 * t211 + (-qJD(3) * t201 * t203 - t229) * t202;
	t197 = t202 * t230 + (t211 - t217) * t200;
	t196 = t207 * t203 + t205 * t229;
	t195 = t207 * t200 - t205 * t230;
	t1 = [t196 * r_i_i_C(1) - t195 * r_i_i_C(2) - t231 * t202 + t206 * t205, t223, (t205 * t213 + t210 * t222) * t204 + (t227 * t202 + t205 * t209) * t201, -t204 * t223 + t217, t199, t197 * r_i_i_C(1) + t198 * r_i_i_C(2); -t198 * r_i_i_C(1) + t197 * r_i_i_C(2) + t206 * t202 + t231 * t205, t199, (t202 * t213 - t210 * t220) * t204 + (t202 * t209 - t227 * t205) * t201, -t202 * t224 - t216, -t223, t195 * r_i_i_C(1) + t196 * r_i_i_C(2); 0, 0, t208 * t204 + (-t210 * t201 - t228) * qJD(3), t221, 0, (t200 * t219 - t203 * t221) * r_i_i_C(2) + (-t200 * t221 - t203 * t219) * r_i_i_C(1);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end