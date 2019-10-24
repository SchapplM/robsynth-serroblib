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
% Datum: 2019-10-24 10:49
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
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
	% StartTime: 2019-10-24 10:49:50
	% EndTime: 2019-10-24 10:49:50
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:49:50
	% EndTime: 2019-10-24 10:49:50
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t6 = cos(qJ(1));
	t5 = sin(qJ(1));
	t1 = [0, 0, 0, 0, 0; (r_i_i_C(1) * t5 + r_i_i_C(2) * t6) * qJD(1), 0, 0, 0, 0; (-r_i_i_C(1) * t6 + r_i_i_C(2) * t5) * qJD(1), 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:49:50
	% EndTime: 2019-10-24 10:49:50
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (22->6), mult. (20->9), div. (0->0), fcn. (10->4), ass. (0->8)
	t15 = qJD(1) + qJD(2);
	t16 = qJ(1) + qJ(2);
	t21 = sin(t16) * t15;
	t20 = cos(t16) * t15;
	t19 = r_i_i_C(1) * t21 + r_i_i_C(2) * t20;
	t18 = pkin(1) * qJD(1);
	t17 = -r_i_i_C(1) * t20 + r_i_i_C(2) * t21;
	t1 = [0, 0, 0, 0, 0; sin(qJ(1)) * t18 + t19, t19, 0, 0, 0; -cos(qJ(1)) * t18 + t17, t17, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:49:50
	% EndTime: 2019-10-24 10:49:50
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (71->13), mult. (58->17), div. (0->0), fcn. (36->6), ass. (0->13)
	t69 = pkin(2) + r_i_i_C(1) * cos(pkin(9));
	t67 = r_i_i_C(2) * sin(pkin(9));
	t58 = qJ(1) + qJ(2);
	t55 = sin(t58);
	t57 = qJD(1) + qJD(2);
	t66 = t57 * t55;
	t56 = cos(t58);
	t65 = t57 * t56;
	t64 = -r_i_i_C(3) - qJ(3);
	t63 = pkin(1) * qJD(1);
	t62 = t65 * t67 + t56 * qJD(3) + (t64 * t55 - t69 * t56) * t57;
	t61 = -t55 * qJD(3) + (-t55 * t67 + t64 * t56) * t57 + t69 * t66;
	t1 = [0, 0, 0, 0, 0; sin(qJ(1)) * t63 + t61, t61, -t66, 0, 0; -cos(qJ(1)) * t63 + t62, t62, t65, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:49:51
	% EndTime: 2019-10-24 10:49:51
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (190->27), mult. (200->44), div. (0->0), fcn. (164->8), ass. (0->25)
	t143 = sin(pkin(9));
	t144 = cos(pkin(9));
	t162 = pkin(3) * t144 + (pkin(7) + r_i_i_C(3)) * t143 + pkin(2);
	t158 = pkin(1) * qJD(1);
	t142 = qJ(1) + qJ(2);
	t139 = sin(t142);
	t141 = qJD(1) + qJD(2);
	t157 = t141 * t139;
	t140 = cos(t142);
	t156 = t141 * t140;
	t145 = sin(qJ(4));
	t155 = t144 * t145;
	t146 = cos(qJ(4));
	t154 = t144 * t146;
	t152 = t139 * t145 + t140 * t154;
	t151 = -t139 * t146 + t140 * t155;
	t150 = t139 * t154 - t140 * t145;
	t149 = t139 * t155 + t140 * t146;
	t130 = -t152 * qJD(4) + t149 * t141;
	t131 = t151 * qJD(4) + t150 * t141;
	t148 = t131 * r_i_i_C(1) - t130 * r_i_i_C(2) - qJ(3) * t156 - t139 * qJD(3) + t162 * t157;
	t132 = t150 * qJD(4) + t151 * t141;
	t133 = -t149 * qJD(4) + t152 * t141;
	t147 = t132 * r_i_i_C(2) - t133 * r_i_i_C(1) + t140 * qJD(3) + (-qJ(3) * t139 - t162 * t140) * t141;
	t1 = [0, 0, 0, (-r_i_i_C(1) * t146 + r_i_i_C(2) * t145) * t143 * qJD(4), 0; sin(qJ(1)) * t158 + t148, t148, -t157, t132 * r_i_i_C(1) + t133 * r_i_i_C(2), 0; -cos(qJ(1)) * t158 + t147, t147, t156, t130 * r_i_i_C(1) + t131 * r_i_i_C(2), 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:49:51
	% EndTime: 2019-10-24 10:49:51
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (390->45), mult. (325->66), div. (0->0), fcn. (262->10), ass. (0->42)
	t196 = cos(qJ(4));
	t217 = pkin(4) * qJD(4);
	t207 = t196 * t217;
	t190 = qJD(1) + qJD(2);
	t193 = sin(pkin(9));
	t214 = t190 * t193;
	t221 = qJD(3) + (-pkin(8) - pkin(7)) * t214 + t207;
	t192 = qJ(1) + qJ(2);
	t186 = sin(t192);
	t189 = qJD(4) + qJD(5);
	t194 = cos(pkin(9));
	t203 = t189 * t194 - t190;
	t220 = t186 * t203;
	t188 = cos(t192);
	t219 = t188 * t203;
	t218 = pkin(1) * qJD(1);
	t216 = t190 * t186;
	t215 = t190 * t188;
	t213 = t190 * t194;
	t195 = sin(qJ(4));
	t212 = t194 * t195;
	t211 = t194 * t196;
	t191 = qJ(4) + qJ(5);
	t185 = sin(t191);
	t187 = cos(t191);
	t204 = -t189 + t213;
	t201 = t204 * t186;
	t171 = t185 * t201 - t187 * t219;
	t172 = t185 * t219 + t187 * t201;
	t210 = t171 * r_i_i_C(1) + t172 * r_i_i_C(2);
	t200 = t204 * t188;
	t173 = t185 * t200 + t187 * t220;
	t174 = -t185 * t220 + t187 * t200;
	t209 = t173 * r_i_i_C(1) + t174 * r_i_i_C(2);
	t208 = r_i_i_C(1) * t187 * t189;
	t205 = -pkin(4) * t195 - qJ(3);
	t202 = t212 * t217;
	t184 = t196 * pkin(4) + pkin(3);
	t199 = -t174 * r_i_i_C(1) + t173 * r_i_i_C(2) + (t205 * t190 + t202) * t186 + ((-r_i_i_C(3) * t193 - t184 * t194 - pkin(2)) * t190 + t221) * t188;
	t198 = pkin(2) * t216 + t172 * r_i_i_C(1) - t171 * r_i_i_C(2) + t188 * t202 + t205 * t215 + (r_i_i_C(3) * t214 + t184 * t213 - t221) * t186;
	t179 = t193 * t189 * t185 * r_i_i_C(2);
	t1 = [0, 0, 0, t179 + (-t207 - t208) * t193, -t193 * t208 + t179; sin(qJ(1)) * t218 + t198, t198, -t216, ((-t186 * t196 + t188 * t212) * t190 + (t186 * t211 - t188 * t195) * qJD(4)) * pkin(4) + t209, t209; -cos(qJ(1)) * t218 + t199, t199, t215, ((t186 * t212 + t188 * t196) * t190 + (-t186 * t195 - t188 * t211) * qJD(4)) * pkin(4) + t210, t210;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end