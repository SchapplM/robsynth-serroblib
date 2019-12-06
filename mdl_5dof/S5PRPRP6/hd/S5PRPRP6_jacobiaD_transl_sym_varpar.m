% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5PRPRP6
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:41
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5PRPRP6_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP6_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP6_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5PRPRP6_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRPRP6_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRPRP6_jacobiaD_transl_sym_varpar: pkin has to be [7x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:41:38
	% EndTime: 2019-12-05 15:41:38
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:41:38
	% EndTime: 2019-12-05 15:41:38
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:41:38
	% EndTime: 2019-12-05 15:41:38
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->2), mult. (16->8), div. (0->0), fcn. (10->4), ass. (0->4)
	t12 = sin(qJ(2));
	t13 = cos(qJ(2));
	t14 = qJD(2) * (-r_i_i_C(1) * t13 + r_i_i_C(2) * t12);
	t1 = [0, cos(pkin(7)) * t14, 0, 0, 0; 0, sin(pkin(7)) * t14, 0, 0, 0; 0, (-r_i_i_C(1) * t12 - r_i_i_C(2) * t13) * qJD(2), 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:41:38
	% EndTime: 2019-12-05 15:41:38
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (12->6), mult. (42->14), div. (0->0), fcn. (30->4), ass. (0->9)
	t97 = -pkin(2) + r_i_i_C(2);
	t96 = r_i_i_C(3) + qJ(3);
	t93 = cos(qJ(2));
	t95 = qJD(2) * t93;
	t92 = sin(qJ(2));
	t94 = qJD(3) * t93 + (-t96 * t92 + t97 * t93) * qJD(2);
	t91 = cos(pkin(7));
	t90 = sin(pkin(7));
	t1 = [0, t94 * t91, t91 * t95, 0, 0; 0, t94 * t90, t90 * t95, 0, 0; 0, t92 * qJD(3) + (t97 * t92 + t96 * t93) * qJD(2), qJD(2) * t92, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:41:38
	% EndTime: 2019-12-05 15:41:38
	% DurationCPUTime: 0.07s
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
	t130 = cos(pkin(7));
	t129 = sin(pkin(7));
	t1 = [0, t135 * t130, t130 * t142, t130 * t137 + ((-t129 * t133 - t130 * t145) * r_i_i_C(1) + (t129 * t131 - t130 * t144) * r_i_i_C(2)) * qJD(4), 0; 0, t135 * t129, t129 * t142, t129 * t137 + ((-t129 * t145 + t130 * t133) * r_i_i_C(1) + (-t129 * t144 - t130 * t131) * r_i_i_C(2)) * qJD(4), 0; 0, t136 * t132 + (t140 * t132 + t138 * t134) * qJD(2), t143, (-t131 * t143 + t133 * t141) * r_i_i_C(2) + (t131 * t141 + t133 * t143) * r_i_i_C(1), 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:41:39
	% EndTime: 2019-12-05 15:41:39
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (78->30), mult. (264->51), div. (0->0), fcn. (214->6), ass. (0->23)
	t181 = sin(qJ(4));
	t183 = cos(qJ(4));
	t196 = r_i_i_C(3) + qJ(5);
	t197 = pkin(4) + r_i_i_C(1);
	t198 = t197 * t181 - t196 * t183 + qJ(3);
	t182 = sin(qJ(2));
	t195 = t181 * t182;
	t194 = t182 * t183;
	t193 = qJD(2) * t182;
	t184 = cos(qJ(2));
	t192 = qJD(2) * t184;
	t191 = qJD(4) * t184;
	t179 = sin(pkin(7));
	t190 = t179 * t192;
	t180 = cos(pkin(7));
	t189 = t180 * t192;
	t188 = t179 * t183 + t180 * t195;
	t187 = -t179 * t195 + t180 * t183;
	t186 = -qJD(5) * t183 + qJD(3) + (-pkin(2) - pkin(6) - r_i_i_C(2)) * qJD(2) + (t196 * t181 + t197 * t183) * qJD(4);
	t185 = t186 * t184 - t198 * t193;
	t177 = t187 * qJD(4) + t183 * t190;
	t175 = t188 * qJD(4) - t183 * t189;
	t1 = [0, t185 * t180, t189, t188 * qJD(5) - t197 * t175 - t196 * (-t181 * t189 + (t179 * t181 - t180 * t194) * qJD(4)), t175; 0, t185 * t179, t190, -t187 * qJD(5) + t197 * t177 + t196 * (t181 * t190 + (t179 * t194 + t180 * t181) * qJD(4)), -t177; 0, t186 * t182 + t198 * t192, t193, (-t196 * t191 + t197 * t193) * t183 + (t196 * t193 + (t197 * qJD(4) - qJD(5)) * t184) * t181, -t181 * t191 - t183 * t193;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end