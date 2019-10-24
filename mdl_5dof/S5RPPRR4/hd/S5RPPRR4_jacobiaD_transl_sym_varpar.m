% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RPPRR4
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:41
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RPPRR4_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR4_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RPPRR4_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPRR4_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:41:09
	% EndTime: 2019-10-24 10:41:09
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:41:09
	% EndTime: 2019-10-24 10:41:09
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
	% StartTime: 2019-10-24 10:41:09
	% EndTime: 2019-10-24 10:41:09
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (11->8), mult. (28->12), div. (0->0), fcn. (18->4), ass. (0->5)
	t37 = -r_i_i_C(3) - qJ(2);
	t36 = r_i_i_C(1) * cos(pkin(8)) - r_i_i_C(2) * sin(pkin(8)) + pkin(1);
	t35 = cos(qJ(1));
	t34 = sin(qJ(1));
	t1 = [0, 0, 0, 0, 0; -t34 * qJD(2) + (t36 * t34 + t37 * t35) * qJD(1), -qJD(1) * t34, 0, 0, 0; t35 * qJD(2) + (t37 * t34 - t36 * t35) * qJD(1), qJD(1) * t35, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:41:09
	% EndTime: 2019-10-24 10:41:09
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (23->16), mult. (62->21), div. (0->0), fcn. (48->6), ass. (0->11)
	t75 = sin(pkin(9));
	t76 = sin(pkin(8));
	t77 = cos(pkin(9));
	t87 = (r_i_i_C(3) + qJ(3)) * t76 + (r_i_i_C(1) * t77 - r_i_i_C(2) * t75 + pkin(2)) * cos(pkin(8)) + pkin(1);
	t79 = sin(qJ(1));
	t85 = qJD(1) * t79;
	t80 = cos(qJ(1));
	t84 = qJD(1) * t80;
	t83 = t76 * qJD(3);
	t81 = -t75 * r_i_i_C(1) - t77 * r_i_i_C(2) - qJ(2);
	t1 = [0, 0, 0, 0, 0; -t80 * t83 - t79 * qJD(2) + (t87 * t79 + t81 * t80) * qJD(1), -t85, -t76 * t84, 0, 0; -t79 * t83 + t80 * qJD(2) + (t81 * t79 - t87 * t80) * qJD(1), t84, -t76 * t85, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:41:09
	% EndTime: 2019-10-24 10:41:09
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (86->31), mult. (140->49), div. (0->0), fcn. (118->8), ass. (0->23)
	t125 = sin(pkin(8));
	t126 = cos(pkin(8));
	t142 = (r_i_i_C(3) + pkin(6) + qJ(3)) * t125 + (cos(pkin(9)) * pkin(3) + pkin(2)) * t126 + pkin(1);
	t128 = sin(qJ(1));
	t140 = t126 * t128;
	t129 = cos(qJ(1));
	t139 = t126 * t129;
	t138 = qJD(1) * t128;
	t137 = qJD(1) * t129;
	t136 = qJD(3) * t125;
	t134 = -sin(pkin(9)) * pkin(3) - qJ(2);
	t123 = pkin(9) + qJ(4);
	t121 = sin(t123);
	t122 = cos(t123);
	t133 = t121 * t128 + t122 * t139;
	t132 = -t121 * t129 + t122 * t140;
	t131 = t121 * t139 - t122 * t128;
	t130 = t121 * t140 + t122 * t129;
	t119 = t133 * qJD(1) - t130 * qJD(4);
	t118 = t131 * qJD(1) + t132 * qJD(4);
	t117 = t132 * qJD(1) + t131 * qJD(4);
	t116 = t130 * qJD(1) - t133 * qJD(4);
	t1 = [0, 0, 0, (-r_i_i_C(1) * t122 + r_i_i_C(2) * t121) * t125 * qJD(4), 0; -t129 * t136 + t117 * r_i_i_C(1) - t116 * r_i_i_C(2) - t128 * qJD(2) + (t142 * t128 + t134 * t129) * qJD(1), -t138, -t125 * t137, t118 * r_i_i_C(1) + t119 * r_i_i_C(2), 0; -t128 * t136 - t119 * r_i_i_C(1) + t118 * r_i_i_C(2) + t129 * qJD(2) + (t134 * t128 - t142 * t129) * qJD(1), t137, -t125 * t138, t116 * r_i_i_C(1) + t117 * r_i_i_C(2), 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:41:09
	% EndTime: 2019-10-24 10:41:09
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (235->48), mult. (239->69), div. (0->0), fcn. (198->10), ass. (0->34)
	t163 = pkin(9) + qJ(4);
	t160 = cos(t163);
	t165 = sin(pkin(8));
	t166 = cos(pkin(8));
	t187 = (r_i_i_C(3) + pkin(7) + pkin(6) + qJ(3)) * t165 + (pkin(4) * t160 + cos(pkin(9)) * pkin(3) + pkin(2)) * t166 + pkin(1);
	t185 = pkin(4) * qJD(4);
	t167 = sin(qJ(1));
	t184 = t166 * t167;
	t168 = cos(qJ(1));
	t183 = t166 * t168;
	t159 = sin(t163);
	t182 = -qJ(2) - pkin(4) * t159 - sin(pkin(9)) * pkin(3);
	t164 = qJD(4) + qJD(5);
	t161 = qJ(5) + t163;
	t157 = sin(t161);
	t158 = cos(t161);
	t171 = t157 * t184 + t158 * t168;
	t174 = t157 * t167 + t158 * t183;
	t150 = t171 * qJD(1) - t164 * t174;
	t172 = t157 * t183 - t158 * t167;
	t173 = -t157 * t168 + t158 * t184;
	t151 = qJD(1) * t173 + t164 * t172;
	t181 = t150 * r_i_i_C(1) + t151 * r_i_i_C(2);
	t152 = qJD(1) * t172 + t164 * t173;
	t153 = qJD(1) * t174 - t171 * t164;
	t180 = t152 * r_i_i_C(1) + t153 * r_i_i_C(2);
	t179 = qJD(1) * t167;
	t178 = qJD(1) * t168;
	t177 = qJD(3) * t165;
	t176 = r_i_i_C(1) * t158 * t164;
	t170 = t159 * t183 - t160 * t167;
	t169 = t159 * t184 + t160 * t168;
	t154 = t165 * t164 * t157 * r_i_i_C(2);
	t1 = [0, 0, 0, t154 + (-t160 * t185 - t176) * t165, -t165 * t176 + t154; -t168 * t177 + t151 * r_i_i_C(1) - t150 * r_i_i_C(2) - t167 * qJD(2) + t170 * t185 + (t187 * t167 + t182 * t168) * qJD(1), -t179, -t165 * t178, ((-t159 * t168 + t160 * t184) * qJD(4) + t170 * qJD(1)) * pkin(4) + t180, t180; -t167 * t177 - t153 * r_i_i_C(1) + t152 * r_i_i_C(2) + t168 * qJD(2) + t169 * t185 + (t182 * t167 - t187 * t168) * qJD(1), t178, -t165 * t179, ((-t159 * t167 - t160 * t183) * qJD(4) + t169 * qJD(1)) * pkin(4) + t181, t181;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end