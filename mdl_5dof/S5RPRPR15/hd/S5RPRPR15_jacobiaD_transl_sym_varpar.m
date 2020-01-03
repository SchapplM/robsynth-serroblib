% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RPRPR15
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RPRPR15_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR15_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR15_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RPRPR15_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRPR15_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR15_jacobiaD_transl_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:37:52
	% EndTime: 2019-12-31 18:37:52
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:37:52
	% EndTime: 2019-12-31 18:37:52
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:37:52
	% EndTime: 2019-12-31 18:37:52
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->6), mult. (20->10), div. (0->0), fcn. (12->2), ass. (0->5)
	t10 = -pkin(1) + r_i_i_C(2);
	t9 = r_i_i_C(3) + qJ(2);
	t8 = cos(qJ(1));
	t7 = sin(qJ(1));
	t1 = [t8 * qJD(2) + (t10 * t8 - t9 * t7) * qJD(1), qJD(1) * t8, 0, 0, 0; t7 * qJD(2) + (t10 * t7 + t9 * t8) * qJD(1), qJD(1) * t7, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:37:52
	% EndTime: 2019-12-31 18:37:53
	% DurationCPUTime: 0.08s
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
	t24 = -pkin(1) - pkin(6) - r_i_i_C(3);
	t22 = r_i_i_C(1) * t17 + r_i_i_C(2) * t19 + qJ(2);
	t21 = qJD(2) + t29;
	t1 = [t21 * t20 + (-t18 * t22 + t20 * t24) * qJD(1), t27, (-t17 * t27 - t19 * t26) * r_i_i_C(2) + (-t17 * t26 + t19 * t27) * r_i_i_C(1), 0, 0; t21 * t18 + (t18 * t24 + t20 * t22) * qJD(1), t28, (-t17 * t28 + t19 * t25) * r_i_i_C(2) + (t17 * t25 + t19 * t28) * r_i_i_C(1), 0, 0; 0, 0, -t29, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:37:53
	% EndTime: 2019-12-31 18:37:53
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (59->23), mult. (186->39), div. (0->0), fcn. (139->6), ass. (0->19)
	t164 = sin(qJ(3));
	t166 = cos(qJ(3));
	t162 = sin(pkin(8));
	t163 = cos(pkin(8));
	t170 = r_i_i_C(1) * t163 - r_i_i_C(2) * t162 + pkin(3);
	t177 = r_i_i_C(3) + qJ(4);
	t181 = -(t177 * t164 + t170 * t166) * qJD(3) + t166 * qJD(4);
	t180 = t170 * qJD(3) - qJD(4);
	t178 = t170 * t164 - t177 * t166 + qJ(2);
	t165 = sin(qJ(1));
	t176 = qJD(1) * t165;
	t167 = cos(qJ(1));
	t175 = qJD(1) * t167;
	t174 = qJD(3) * t165;
	t173 = qJD(3) * t167;
	t171 = qJD(1) * t177;
	t169 = qJD(1) * t170;
	t168 = qJD(2) + (-t162 * r_i_i_C(1) - t163 * r_i_i_C(2) - pkin(1) - pkin(6)) * qJD(1) - t181;
	t1 = [t168 * t167 - t178 * t176, t175, (t167 * t169 + t177 * t174) * t166 + (-t180 * t165 + t167 * t171) * t164, t164 * t174 - t166 * t175, 0; t168 * t165 + t178 * t175, t176, (t165 * t169 - t177 * t173) * t166 + (t165 * t171 + t180 * t167) * t164, -t164 * t173 - t166 * t176, 0; 0, 0, t181, qJD(3) * t166, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:37:53
	% EndTime: 2019-12-31 18:37:53
	% DurationCPUTime: 0.21s
	% Computational Cost: add. (174->44), mult. (311->71), div. (0->0), fcn. (248->8), ass. (0->37)
	t209 = sin(qJ(3));
	t203 = cos(pkin(8)) * pkin(4) + pkin(3);
	t206 = pkin(8) + qJ(5);
	t204 = sin(t206);
	t205 = cos(t206);
	t216 = r_i_i_C(1) * t205 - r_i_i_C(2) * t204 + t203;
	t211 = cos(qJ(3));
	t236 = r_i_i_C(3) + pkin(7) + qJ(4);
	t239 = t236 * t211;
	t245 = -(-t216 * t209 + t239) * qJD(3) - qJD(4) * t209;
	t224 = t236 * t209;
	t244 = t216 * t211 + t224;
	t243 = -t203 * t209 - qJ(2) + t239;
	t210 = sin(qJ(1));
	t220 = qJD(1) * t209 + qJD(5);
	t212 = cos(qJ(1));
	t231 = qJD(3) * t212;
	t238 = t220 * t210 - t211 * t231;
	t232 = qJD(3) * t211;
	t237 = t210 * t232 + t220 * t212;
	t235 = qJD(1) * t210;
	t234 = qJD(1) * t212;
	t233 = qJD(3) * t209;
	t229 = qJD(5) * t211;
	t228 = t211 * qJD(4);
	t221 = -qJD(5) * t209 - qJD(1);
	t219 = r_i_i_C(1) * t204 + r_i_i_C(2) * t205;
	t218 = t221 * t210;
	t217 = t221 * t212;
	t215 = qJD(5) * t219;
	t214 = qJD(1) * t244;
	t213 = -t228 + qJD(2) + (t203 * t211 + t224) * qJD(3) + (-pkin(4) * sin(pkin(8)) - pkin(1) - pkin(6)) * qJD(1);
	t202 = t204 * t218 + t237 * t205;
	t201 = -t237 * t204 + t205 * t218;
	t200 = t204 * t217 - t238 * t205;
	t199 = t238 * t204 + t205 * t217;
	t1 = [t200 * r_i_i_C(1) + t199 * r_i_i_C(2) + t213 * t212 + t243 * t235, t234, t212 * t214 + (-t219 * t229 - t245) * t210, t210 * t233 - t211 * t234, r_i_i_C(1) * t201 - r_i_i_C(2) * t202; t202 * r_i_i_C(1) + t201 * r_i_i_C(2) + t213 * t210 - t243 * t234, t235, t210 * t214 + (t211 * t215 + t245) * t212, -t209 * t231 - t211 * t235, -r_i_i_C(1) * t199 + r_i_i_C(2) * t200; 0, 0, -t244 * qJD(3) + t209 * t215 + t228, t232, (t204 * t229 + t205 * t233) * r_i_i_C(2) + (t204 * t233 - t205 * t229) * r_i_i_C(1);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end