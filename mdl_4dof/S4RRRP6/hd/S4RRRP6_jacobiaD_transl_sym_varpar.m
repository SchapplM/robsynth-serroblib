% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S4RRRP6
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% JaD_transl [3x4]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 14:02
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S4RRRP6_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),uint8(0),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP6_jacobiaD_transl_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP6_jacobiaD_transl_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S4RRRP6_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4RRRP6_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP6_jacobiaD_transl_sym_varpar: pkin has to be [6x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 14:02:09
	% EndTime: 2019-12-29 14:02:09
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 14:02:14
	% EndTime: 2019-12-29 14:02:14
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0; 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 14:02:14
	% EndTime: 2019-12-29 14:02:14
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (19->15), mult. (64->29), div. (0->0), fcn. (40->4), ass. (0->13)
	t28 = pkin(5) + r_i_i_C(3);
	t18 = sin(qJ(1));
	t27 = qJD(1) * t18;
	t20 = cos(qJ(1));
	t26 = qJD(1) * t20;
	t25 = qJD(2) * t18;
	t24 = qJD(2) * t20;
	t17 = sin(qJ(2));
	t19 = cos(qJ(2));
	t23 = r_i_i_C(1) * t17 + r_i_i_C(2) * t19;
	t22 = -r_i_i_C(1) * t19 + r_i_i_C(2) * t17 - pkin(1);
	t21 = t23 * qJD(2);
	t1 = [t23 * t25 + (-t18 * t28 + t20 * t22) * qJD(1), (t17 * t24 + t19 * t27) * r_i_i_C(2) + (t17 * t27 - t19 * t24) * r_i_i_C(1), 0, 0; -t20 * t21 + (t18 * t22 + t20 * t28) * qJD(1), (t17 * t25 - t19 * t26) * r_i_i_C(2) + (-t17 * t26 - t19 * t25) * r_i_i_C(1), 0, 0; 0, -t21, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 14:02:16
	% EndTime: 2019-12-29 14:02:16
	% DurationCPUTime: 0.29s
	% Computational Cost: add. (83->38), mult. (270->73), div. (0->0), fcn. (211->6), ass. (0->31)
	t199 = sin(qJ(2));
	t202 = cos(qJ(2));
	t223 = pkin(6) + r_i_i_C(3);
	t213 = t223 * t202;
	t224 = -pkin(2) * t199 + t213;
	t201 = cos(qJ(3));
	t203 = cos(qJ(1));
	t221 = t201 * t203;
	t200 = sin(qJ(1));
	t220 = qJD(1) * t200;
	t219 = qJD(1) * t203;
	t218 = qJD(2) * t200;
	t217 = qJD(2) * t202;
	t216 = qJD(2) * t203;
	t215 = qJD(3) * t199;
	t214 = qJD(3) * t202;
	t212 = -qJD(1) + t214;
	t211 = qJD(1) * t202 - qJD(3);
	t198 = sin(qJ(3));
	t210 = r_i_i_C(1) * t198 + r_i_i_C(2) * t201;
	t209 = r_i_i_C(1) * t201 - r_i_i_C(2) * t198 + pkin(2);
	t208 = t212 * t198;
	t207 = -pkin(2) * t202 - t223 * t199 - pkin(1);
	t206 = qJD(2) * t209;
	t205 = t199 * t216 + t211 * t200;
	t204 = -t223 * qJD(2) + t210 * qJD(3);
	t197 = -t211 * t221 + (qJD(2) * t199 * t201 + t208) * t200;
	t196 = t212 * t201 * t200 + (-t199 * t218 + t211 * t203) * t198;
	t195 = t205 * t201 + t203 * t208;
	t194 = t205 * t198 - t212 * t221;
	t1 = [t197 * r_i_i_C(1) + t196 * r_i_i_C(2) - t224 * t218 + (-pkin(5) * t200 + t207 * t203) * qJD(1), (-t203 * t206 - t223 * t220) * t202 + (t204 * t203 + t209 * t220) * t199, t194 * r_i_i_C(1) + t195 * r_i_i_C(2), 0; -t195 * r_i_i_C(1) + t194 * r_i_i_C(2) + t224 * t216 + (pkin(5) * t203 + t207 * t200) * qJD(1), (-t200 * t206 + t223 * t219) * t202 + (t204 * t200 - t209 * t219) * t199, -t196 * r_i_i_C(1) + t197 * r_i_i_C(2), 0; 0, -t210 * t214 + (-t209 * t199 + t213) * qJD(2), (t198 * t215 - t201 * t217) * r_i_i_C(2) + (-t198 * t217 - t201 * t215) * r_i_i_C(1), 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 14:02:16
	% EndTime: 2019-12-29 14:02:16
	% DurationCPUTime: 0.40s
	% Computational Cost: add. (129->42), mult. (370->69), div. (0->0), fcn. (290->6), ass. (0->38)
	t237 = pkin(3) + r_i_i_C(1);
	t203 = sin(qJ(3));
	t204 = sin(qJ(2));
	t205 = sin(qJ(1));
	t208 = cos(qJ(1));
	t207 = cos(qJ(2));
	t219 = qJD(1) * t207 - qJD(3);
	t215 = t219 * t208;
	t206 = cos(qJ(3));
	t220 = qJD(3) * t207 - qJD(1);
	t217 = t220 * t206;
	t228 = qJD(2) * t205;
	t199 = (-t204 * t228 + t215) * t203 + t205 * t217;
	t234 = t206 * pkin(3);
	t201 = pkin(2) + t234;
	t216 = t220 * t203;
	t225 = t204 * qJD(4);
	t233 = r_i_i_C(3) + qJ(4) + pkin(6);
	t239 = t233 * t207;
	t246 = (-t201 * t204 + t239) * qJD(2) - pkin(3) * t216 + pkin(5) * qJD(1) + t225;
	t236 = r_i_i_C(2) * t203;
	t214 = r_i_i_C(1) * t206 + t201 - t236;
	t211 = -t214 * t204 + t239;
	t238 = r_i_i_C(2) * t206 + t237 * t203;
	t231 = qJD(1) * t205;
	t230 = qJD(1) * t208;
	t229 = qJD(2) * t204;
	t227 = qJD(2) * t208;
	t226 = qJD(3) * t204;
	t223 = t233 * t204;
	t213 = t238 * t207;
	t212 = t204 * t227 + t219 * t205;
	t210 = qJD(3) * t234 + (-t201 * t207 - pkin(1) - t223) * qJD(1);
	t197 = t212 * t203 - t208 * t217;
	t209 = qJD(4) * t207 + t238 * t226 + (-t214 * t207 - t223) * qJD(2);
	t200 = -t206 * t215 + (t206 * t229 + t216) * t205;
	t198 = t212 * t206 + t208 * t216;
	t1 = [t200 * r_i_i_C(1) + t199 * r_i_i_C(2) - t246 * t205 + t210 * t208, t209 * t208 - t211 * t231, t198 * r_i_i_C(2) + t237 * t197, -t204 * t231 + t207 * t227; -t198 * r_i_i_C(1) + t197 * r_i_i_C(2) + t210 * t205 + t246 * t208, t209 * t205 + t211 * t230, t200 * r_i_i_C(2) - t237 * t199, t204 * t230 + t207 * t228; 0, t211 * qJD(2) - qJD(3) * t213 + t225, (-t237 * t206 + t236) * t226 - qJD(2) * t213, t229;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,4);
end