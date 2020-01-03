% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RPPPR6
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RPPPR6_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR6_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR6_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RPPPR6_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPPR6_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR6_jacobiaD_transl_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:48:16
	% EndTime: 2019-12-31 17:48:16
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:48:16
	% EndTime: 2019-12-31 17:48:16
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
	% StartTime: 2019-12-31 17:48:16
	% EndTime: 2019-12-31 17:48:16
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (10->7), mult. (28->12), div. (0->0), fcn. (18->4), ass. (0->5)
	t16 = r_i_i_C(3) + qJ(2);
	t15 = -r_i_i_C(1) * cos(pkin(7)) + r_i_i_C(2) * sin(pkin(7)) - pkin(1);
	t14 = cos(qJ(1));
	t13 = sin(qJ(1));
	t1 = [t14 * qJD(2) + (-t16 * t13 + t15 * t14) * qJD(1), qJD(1) * t14, 0, 0, 0; t13 * qJD(2) + (t15 * t13 + t16 * t14) * qJD(1), qJD(1) * t13, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:48:16
	% EndTime: 2019-12-31 17:48:16
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (17->12), mult. (48->17), div. (0->0), fcn. (34->4), ass. (0->9)
	t98 = r_i_i_C(1) + qJ(2);
	t92 = sin(qJ(1));
	t97 = qJD(1) * t92;
	t93 = cos(qJ(1));
	t96 = qJD(1) * t93;
	t90 = sin(pkin(7));
	t95 = qJD(3) * t90;
	t94 = -pkin(1) + (-pkin(2) + r_i_i_C(2)) * cos(pkin(7)) + (-r_i_i_C(3) - qJ(3)) * t90;
	t1 = [-t92 * t95 + t93 * qJD(2) + (-t98 * t92 + t94 * t93) * qJD(1), t96, -t90 * t97, 0, 0; t93 * t95 + t92 * qJD(2) + (t94 * t92 + t98 * t93) * qJD(1), t97, t90 * t96, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:48:16
	% EndTime: 2019-12-31 17:48:16
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (28->18), mult. (80->24), div. (0->0), fcn. (62->6), ass. (0->12)
	t117 = sin(qJ(1));
	t123 = qJD(1) * t117;
	t118 = cos(qJ(1));
	t122 = qJD(1) * t118;
	t114 = sin(pkin(7));
	t116 = cos(pkin(7));
	t121 = qJD(3) * t114 + qJD(4) * t116;
	t113 = sin(pkin(8));
	t115 = cos(pkin(8));
	t120 = r_i_i_C(1) * t115 - r_i_i_C(2) * t113 + pkin(3) + qJ(2);
	t119 = -pkin(1) + (-pkin(2) - r_i_i_C(3) - qJ(4)) * t116 + (-r_i_i_C(1) * t113 - r_i_i_C(2) * t115 - qJ(3)) * t114;
	t1 = [t118 * qJD(2) - t121 * t117 + (-t117 * t120 + t118 * t119) * qJD(1), t122, -t114 * t123, -t116 * t123, 0; t117 * qJD(2) + t121 * t118 + (t117 * t119 + t118 * t120) * qJD(1), t123, t114 * t122, t116 * t122, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:48:16
	% EndTime: 2019-12-31 17:48:17
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (79->48), mult. (242->79), div. (0->0), fcn. (222->8), ass. (0->33)
	t227 = r_i_i_C(3) + pkin(6);
	t226 = -pkin(2) - qJ(4);
	t225 = pkin(3) + qJ(2);
	t202 = sin(pkin(8));
	t207 = sin(qJ(1));
	t224 = t202 * t207;
	t209 = cos(qJ(1));
	t223 = t202 * t209;
	t205 = cos(pkin(7));
	t206 = sin(qJ(5));
	t222 = t205 * t206;
	t208 = cos(qJ(5));
	t221 = t205 * t208;
	t220 = t205 * t209;
	t204 = cos(pkin(8));
	t219 = t207 * t204;
	t218 = t209 * t204;
	t217 = qJD(1) * t207;
	t216 = qJD(1) * t209;
	t215 = t205 * t217;
	t214 = t205 * t216;
	t203 = sin(pkin(7));
	t213 = -qJ(3) * t203 - pkin(1);
	t212 = t203 * qJD(3) + t205 * qJD(4);
	t200 = -t203 * t224 + t218;
	t211 = -t200 * t206 - t207 * t221;
	t210 = -t200 * t208 + t207 * t222;
	t199 = t203 * t223 + t219;
	t198 = t200 * qJD(1);
	t196 = t199 * qJD(1);
	t194 = -t206 * t215 + t198 * t208 + (-t199 * t206 + t208 * t220) * qJD(5);
	t193 = -t208 * t215 - t198 * t206 + (-t199 * t208 - t206 * t220) * qJD(5);
	t1 = [t209 * qJD(2) - t212 * t207 + (-t208 * r_i_i_C(1) + t206 * r_i_i_C(2) - pkin(4)) * t196 + (t211 * r_i_i_C(1) + t210 * r_i_i_C(2)) * qJD(5) + (-t227 * (-t203 * t218 + t224) - t225 * t207 + ((-t206 * r_i_i_C(1) - t208 * r_i_i_C(2) + t226) * t205 + t213) * t209) * qJD(1), t216, -t203 * t217, -t215, t193 * r_i_i_C(1) - t194 * r_i_i_C(2); t198 * pkin(4) + t194 * r_i_i_C(1) + t193 * r_i_i_C(2) + t207 * qJD(2) + t212 * t209 + (t227 * (t203 * t219 + t223) + t225 * t209 + (t226 * t205 + t213) * t207) * qJD(1), t217, t203 * t216, t214, (-t196 * t206 + t208 * t214) * r_i_i_C(1) + (-t196 * t208 - t206 * t214) * r_i_i_C(2) + (-t210 * r_i_i_C(1) + t211 * r_i_i_C(2)) * qJD(5); 0, 0, 0, 0, ((t202 * t221 - t203 * t206) * r_i_i_C(1) + (-t202 * t222 - t203 * t208) * r_i_i_C(2)) * qJD(5);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end