% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRPRR6
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RRPRR6_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR6_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR6_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRPRR6_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPRR6_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR6_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:06:32
	% EndTime: 2020-01-03 12:06:32
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:06:32
	% EndTime: 2020-01-03 12:06:32
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t6 = cos(qJ(1));
	t5 = sin(qJ(1));
	t1 = [0, 0, 0, 0, 0; (-r_i_i_C(1) * t5 - r_i_i_C(2) * t6) * qJD(1), 0, 0, 0, 0; (r_i_i_C(1) * t6 - r_i_i_C(2) * t5) * qJD(1), 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:06:32
	% EndTime: 2020-01-03 12:06:32
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (22->6), mult. (20->9), div. (0->0), fcn. (10->4), ass. (0->8)
	t19 = pkin(1) * qJD(1);
	t16 = qJ(1) + qJ(2);
	t13 = sin(t16);
	t14 = cos(t16);
	t15 = qJD(1) + qJD(2);
	t18 = (r_i_i_C(1) * t14 - r_i_i_C(2) * t13) * t15;
	t17 = (-r_i_i_C(1) * t13 - r_i_i_C(2) * t14) * t15;
	t1 = [0, 0, 0, 0, 0; -sin(qJ(1)) * t19 + t17, t17, 0, 0, 0; cos(qJ(1)) * t19 + t18, t18, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:06:32
	% EndTime: 2020-01-03 12:06:32
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (71->13), mult. (58->17), div. (0->0), fcn. (36->6), ass. (0->13)
	t76 = pkin(2) + r_i_i_C(1) * cos(pkin(9));
	t63 = qJD(1) + qJD(2);
	t75 = qJD(3) + r_i_i_C(2) * t63 * sin(pkin(9));
	t64 = qJ(1) + qJ(2);
	t61 = sin(t64);
	t73 = t63 * t61;
	t62 = cos(t64);
	t72 = t63 * t62;
	t71 = pkin(1) * qJD(1);
	t70 = qJ(3) * t63;
	t68 = r_i_i_C(3) * t72 + t75 * t61 + t62 * t70 - t76 * t73;
	t67 = r_i_i_C(3) * t73 + t61 * t70 - t75 * t62 + t76 * t72;
	t1 = [0, 0, 0, 0, 0; -sin(qJ(1)) * t71 + t68, t68, t73, 0, 0; cos(qJ(1)) * t71 + t67, t67, -t72, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:06:32
	% EndTime: 2020-01-03 12:06:32
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (190->27), mult. (200->43), div. (0->0), fcn. (164->8), ass. (0->25)
	t160 = sin(pkin(9));
	t161 = cos(pkin(9));
	t179 = pkin(3) * t161 + (pkin(7) + r_i_i_C(3)) * t160 + pkin(2);
	t175 = pkin(1) * qJD(1);
	t159 = qJ(1) + qJ(2);
	t156 = sin(t159);
	t158 = qJD(1) + qJD(2);
	t174 = t158 * t156;
	t157 = cos(t159);
	t173 = t158 * t157;
	t162 = sin(qJ(4));
	t172 = t161 * t162;
	t163 = cos(qJ(4));
	t171 = t161 * t163;
	t169 = t156 * t162 + t157 * t171;
	t168 = t156 * t163 - t157 * t172;
	t167 = -t156 * t171 + t157 * t162;
	t166 = t156 * t172 + t157 * t163;
	t147 = qJD(4) * t167 + t158 * t168;
	t148 = -qJD(4) * t166 + t158 * t169;
	t165 = t148 * r_i_i_C(1) + t147 * r_i_i_C(2) + qJ(3) * t174 - t157 * qJD(3) + t173 * t179;
	t145 = -qJD(4) * t169 + t158 * t166;
	t146 = qJD(4) * t168 + t158 * t167;
	t164 = t146 * r_i_i_C(1) + t145 * r_i_i_C(2) + qJ(3) * t173 + t156 * qJD(3) - t174 * t179;
	t1 = [0, 0, 0, (-r_i_i_C(1) * t163 + r_i_i_C(2) * t162) * t160 * qJD(4), 0; -sin(qJ(1)) * t175 + t164, t164, t174, r_i_i_C(1) * t147 - r_i_i_C(2) * t148, 0; cos(qJ(1)) * t175 + t165, t165, -t173, -r_i_i_C(1) * t145 + r_i_i_C(2) * t146, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:06:32
	% EndTime: 2020-01-03 12:06:32
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (390->45), mult. (325->67), div. (0->0), fcn. (262->10), ass. (0->41)
	t190 = qJD(1) + qJD(2);
	t195 = sin(qJ(4));
	t194 = cos(pkin(9));
	t212 = t194 * t195;
	t217 = pkin(4) * qJD(4);
	t222 = pkin(4) * t190 * t195 - t212 * t217;
	t196 = cos(qJ(4));
	t206 = t196 * t217;
	t193 = sin(pkin(9));
	t214 = t190 * t193;
	t221 = qJD(3) + (-pkin(8) - pkin(7)) * t214 + t206;
	t192 = qJ(1) + qJ(2);
	t186 = sin(t192);
	t189 = qJD(4) + qJD(5);
	t213 = t190 * t194;
	t204 = -t189 + t213;
	t220 = t186 * t204;
	t188 = cos(t192);
	t219 = t188 * t204;
	t218 = pkin(1) * qJD(1);
	t216 = t190 * t186;
	t215 = t190 * t188;
	t211 = t194 * t196;
	t191 = qJ(4) + qJ(5);
	t185 = sin(t191);
	t187 = cos(t191);
	t203 = -t189 * t194 + t190;
	t200 = t203 * t188;
	t169 = t185 * t220 + t187 * t200;
	t170 = t185 * t200 - t187 * t220;
	t210 = -t169 * r_i_i_C(1) + t170 * r_i_i_C(2);
	t201 = t203 * t186;
	t171 = -t185 * t219 + t187 * t201;
	t172 = t185 * t201 + t187 * t219;
	t209 = t171 * r_i_i_C(1) - t172 * r_i_i_C(2);
	t207 = r_i_i_C(1) * t187 * t189;
	t184 = t196 * pkin(4) + pkin(3);
	t199 = t169 * r_i_i_C(2) + t170 * r_i_i_C(1) + qJ(3) * t215 + (-r_i_i_C(3) * t193 - t184 * t194 - pkin(2)) * t216 + t222 * t188 + t221 * t186;
	t198 = pkin(2) * t215 + t172 * r_i_i_C(1) + t171 * r_i_i_C(2) + qJ(3) * t216 + t222 * t186 + (r_i_i_C(3) * t214 + t184 * t213 - t221) * t188;
	t175 = t193 * t189 * t185 * r_i_i_C(2);
	t1 = [0, 0, 0, t175 + (-t206 - t207) * t193, -t193 * t207 + t175; -sin(qJ(1)) * t218 + t199, t199, t216, ((t186 * t196 - t188 * t212) * t190 + (-t186 * t211 + t188 * t195) * qJD(4)) * pkin(4) + t209, t209; cos(qJ(1)) * t218 + t198, t198, -t215, ((-t186 * t212 - t188 * t196) * t190 + (t186 * t195 + t188 * t211) * qJD(4)) * pkin(4) + t210, t210;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end