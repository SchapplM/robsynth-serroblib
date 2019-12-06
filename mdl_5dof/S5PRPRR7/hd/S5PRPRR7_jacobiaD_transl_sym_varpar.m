% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5PRPRR7
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
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5PRPRR7_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR7_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR7_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5PRPRR7_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRPRR7_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRR7_jacobiaD_transl_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:01:10
	% EndTime: 2019-12-05 16:01:10
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:01:10
	% EndTime: 2019-12-05 16:01:10
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:01:10
	% EndTime: 2019-12-05 16:01:11
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->2), mult. (16->8), div. (0->0), fcn. (10->4), ass. (0->4)
	t12 = sin(qJ(2));
	t13 = cos(qJ(2));
	t14 = qJD(2) * (-r_i_i_C(1) * t13 + r_i_i_C(2) * t12);
	t1 = [0, cos(pkin(8)) * t14, 0, 0, 0; 0, sin(pkin(8)) * t14, 0, 0, 0; 0, (-r_i_i_C(1) * t12 - r_i_i_C(2) * t13) * qJD(2), 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:01:11
	% EndTime: 2019-12-05 16:01:11
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (12->6), mult. (42->14), div. (0->0), fcn. (30->4), ass. (0->9)
	t97 = -pkin(2) + r_i_i_C(2);
	t96 = r_i_i_C(3) + qJ(3);
	t93 = cos(qJ(2));
	t95 = qJD(2) * t93;
	t92 = sin(qJ(2));
	t94 = qJD(3) * t93 + (-t96 * t92 + t97 * t93) * qJD(2);
	t91 = cos(pkin(8));
	t90 = sin(pkin(8));
	t1 = [0, t94 * t91, t91 * t95, 0, 0; 0, t94 * t90, t90 * t95, 0, 0; 0, t92 * qJD(3) + (t97 * t92 + t96 * t93) * qJD(2), qJD(2) * t92, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:01:11
	% EndTime: 2019-12-05 16:01:12
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (37->20), mult. (134->45), div. (0->0), fcn. (102->6), ass. (0->18)
	t131 = sin(qJ(4));
	t132 = sin(qJ(2));
	t145 = t131 * t132;
	t133 = cos(qJ(4));
	t144 = t132 * t133;
	t143 = qJD(2) * t132;
	t134 = cos(qJ(2));
	t142 = qJD(2) * t134;
	t141 = qJD(4) * t134;
	t140 = -pkin(2) - pkin(6) - r_i_i_C(3);
	t139 = r_i_i_C(1) * t133 - r_i_i_C(2) * t131;
	t138 = r_i_i_C(1) * t131 + r_i_i_C(2) * t133 + qJ(3);
	t137 = t139 * t142;
	t136 = t139 * qJD(4) + qJD(3);
	t135 = t136 * t134 + (-t138 * t132 + t140 * t134) * qJD(2);
	t130 = cos(pkin(8));
	t129 = sin(pkin(8));
	t1 = [0, t135 * t130, t130 * t142, t130 * t137 + ((-t129 * t133 - t130 * t145) * r_i_i_C(1) + (t129 * t131 - t130 * t144) * r_i_i_C(2)) * qJD(4), 0; 0, t135 * t129, t129 * t142, t129 * t137 + ((-t129 * t145 + t130 * t133) * r_i_i_C(1) + (-t129 * t144 - t130 * t131) * r_i_i_C(2)) * qJD(4), 0; 0, t136 * t132 + (t140 * t132 + t138 * t134) * qJD(2), t143, (-t131 * t143 + t133 * t141) * r_i_i_C(2) + (t131 * t141 + t133 * t143) * r_i_i_C(1), 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:01:12
	% EndTime: 2019-12-05 16:01:12
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (137->30), mult. (232->57), div. (0->0), fcn. (178->8), ass. (0->30)
	t169 = qJ(4) + qJ(5);
	t166 = sin(t169);
	t167 = cos(t169);
	t200 = r_i_i_C(1) * t166 + r_i_i_C(2) * t167;
	t199 = r_i_i_C(1) * t167 - r_i_i_C(2) * t166;
	t168 = qJD(4) + qJD(5);
	t173 = sin(qJ(2));
	t194 = t168 * t173;
	t172 = sin(qJ(4));
	t192 = t172 * t173;
	t170 = sin(pkin(8));
	t171 = cos(pkin(8));
	t175 = cos(qJ(2));
	t188 = qJD(2) * t175;
	t183 = t171 * t188;
	t182 = -t168 * t170 + t183;
	t185 = t171 * t194;
	t191 = (-t166 * t185 + t182 * t167) * r_i_i_C(1) + (-t182 * t166 - t167 * t185) * r_i_i_C(2);
	t184 = t170 * t188;
	t181 = t168 * t171 + t184;
	t186 = t170 * t194;
	t190 = (-t166 * t186 + t181 * t167) * r_i_i_C(1) + (-t181 * t166 - t167 * t186) * r_i_i_C(2);
	t189 = qJD(2) * t173;
	t187 = -pkin(2) - r_i_i_C(3) - pkin(7) - pkin(6);
	t180 = t200 * t168 * t175 + t199 * t189;
	t179 = pkin(4) * t172 + qJ(3) + t200;
	t174 = cos(qJ(4));
	t178 = pkin(4) * qJD(4) * t174 + t199 * t168 + qJD(3);
	t177 = t178 * t175 + (-t179 * t173 + t187 * t175) * qJD(2);
	t1 = [0, t177 * t171, t183, (t174 * t183 + (-t170 * t174 - t171 * t192) * qJD(4)) * pkin(4) + t191, t191; 0, t177 * t170, t184, (t174 * t184 + (-t170 * t192 + t171 * t174) * qJD(4)) * pkin(4) + t190, t190; 0, t178 * t173 + (t187 * t173 + t179 * t175) * qJD(2), t189, (qJD(4) * t172 * t175 + t174 * t189) * pkin(4) + t180, t180;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end