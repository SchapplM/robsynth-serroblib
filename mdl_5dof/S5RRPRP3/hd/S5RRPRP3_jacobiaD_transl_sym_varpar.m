% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRPRP3
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RRPRP3_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP3_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP3_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRPRP3_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPRP3_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP3_jacobiaD_transl_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:51:49
	% EndTime: 2019-12-31 19:51:49
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:51:49
	% EndTime: 2019-12-31 19:51:49
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
	% StartTime: 2019-12-31 19:51:49
	% EndTime: 2019-12-31 19:51:49
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (22->6), mult. (20->9), div. (0->0), fcn. (10->4), ass. (0->8)
	t43 = pkin(1) * qJD(1);
	t40 = qJ(1) + qJ(2);
	t37 = sin(t40);
	t38 = cos(t40);
	t39 = qJD(1) + qJD(2);
	t42 = (-r_i_i_C(1) * t38 + r_i_i_C(2) * t37) * t39;
	t41 = (-r_i_i_C(1) * t37 - r_i_i_C(2) * t38) * t39;
	t1 = [-cos(qJ(1)) * t43 + t42, t42, 0, 0, 0; -sin(qJ(1)) * t43 + t41, t41, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:51:48
	% EndTime: 2019-12-31 19:51:49
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (70->11), mult. (58->14), div. (0->0), fcn. (36->6), ass. (0->13)
	t37 = qJ(1) + qJ(2);
	t34 = sin(t37);
	t36 = qJD(1) + qJD(2);
	t46 = t36 * t34;
	t48 = qJ(3) + r_i_i_C(3);
	t47 = qJD(3) + r_i_i_C(2) * t36 * sin(pkin(8));
	t35 = cos(t37);
	t45 = t36 * t35;
	t44 = pkin(1) * qJD(1);
	t42 = -r_i_i_C(1) * cos(pkin(8)) - pkin(2);
	t41 = t47 * t34 + t42 * t46 + t48 * t45;
	t40 = -t48 * t46 + (t42 * t36 + t47) * t35;
	t1 = [-cos(qJ(1)) * t44 + t40, t40, t45, 0, 0; -sin(qJ(1)) * t44 + t41, t41, t46, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:51:49
	% EndTime: 2019-12-31 19:51:49
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (139->25), mult. (114->37), div. (0->0), fcn. (74->7), ass. (0->21)
	t48 = pkin(8) + qJ(4);
	t44 = sin(t48);
	t45 = cos(t48);
	t50 = qJ(1) + qJ(2);
	t46 = sin(t50);
	t61 = qJD(4) * t46;
	t47 = cos(t50);
	t49 = qJD(1) + qJD(2);
	t64 = t49 * t47;
	t66 = t44 * t64 + t45 * t61;
	t65 = t49 * t46;
	t63 = t49 * (-pkin(7) - qJ(3));
	t62 = pkin(1) * qJD(1);
	t60 = qJD(4) * t47;
	t59 = t44 * t65;
	t57 = t44 * t61;
	t55 = -r_i_i_C(1) * t45 - cos(pkin(8)) * pkin(3) - pkin(2);
	t54 = (-r_i_i_C(1) * t44 - r_i_i_C(2) * t45) * qJD(4);
	t53 = r_i_i_C(1) * t57 + t46 * t63 + t47 * qJD(3) + (-r_i_i_C(3) * t46 + t55 * t47) * t49 + t66 * r_i_i_C(2);
	t52 = r_i_i_C(2) * t59 + r_i_i_C(3) * t64 + t46 * qJD(3) + t55 * t65 + (t54 - t63) * t47;
	t1 = [-cos(qJ(1)) * t62 + t53, t53, t64, (t44 * t60 + t45 * t65) * r_i_i_C(2) + (-t45 * t60 + t59) * r_i_i_C(1), 0; -sin(qJ(1)) * t62 + t52, t52, t65, (-t45 * t64 + t57) * r_i_i_C(2) - t66 * r_i_i_C(1), 0; 0, 0, 0, t54, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:51:49
	% EndTime: 2019-12-31 19:51:49
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (260->34), mult. (212->47), div. (0->0), fcn. (146->7), ass. (0->24)
	t178 = pkin(8) + qJ(4);
	t174 = sin(t178);
	t196 = pkin(4) + r_i_i_C(1);
	t197 = t196 * t174;
	t195 = r_i_i_C(3) + qJ(5);
	t194 = pkin(1) * qJD(1);
	t180 = qJ(1) + qJ(2);
	t176 = sin(t180);
	t179 = qJD(1) + qJD(2);
	t193 = t179 * t176;
	t177 = cos(t180);
	t192 = t179 * t177;
	t191 = qJD(4) * t176;
	t190 = t174 * qJD(5);
	t189 = t196 * qJD(4);
	t175 = cos(t178);
	t187 = qJD(4) * t175 * t177;
	t186 = qJD(4) * t195;
	t185 = qJD(5) - t189;
	t184 = -t195 * t174 - t196 * t175 - cos(pkin(8)) * pkin(3) - pkin(2);
	t181 = -pkin(7) - qJ(3);
	t183 = r_i_i_C(2) * t192 + t195 * t187 + (t184 * t179 + qJD(3)) * t176 + (-t174 * t189 - t179 * t181 + t190) * t177;
	t182 = t181 * t193 + t177 * qJD(3) + (-t175 * t186 - t190) * t176 + (-r_i_i_C(2) * t176 + t184 * t177) * t179 + t191 * t197;
	t1 = [-cos(qJ(1)) * t194 + t182, t182, t192, (-t177 * t186 + t196 * t193) * t174 + (t185 * t177 - t195 * t193) * t175, -t174 * t193 + t187; -sin(qJ(1)) * t194 + t183, t183, t193, (-t176 * t186 - t196 * t192) * t174 + (t185 * t176 + t195 * t192) * t175, t174 * t192 + t175 * t191; 0, 0, 0, t190 + (t195 * t175 - t197) * qJD(4), qJD(4) * t174;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end