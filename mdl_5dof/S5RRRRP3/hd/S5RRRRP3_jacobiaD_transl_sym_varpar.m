% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRRRP3
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RRRRP3_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP3_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP3_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRRRP3_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRP3_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP3_jacobiaD_transl_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:49:51
	% EndTime: 2019-12-31 21:49:51
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:49:51
	% EndTime: 2019-12-31 21:49:51
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
	% StartTime: 2019-12-31 21:49:51
	% EndTime: 2019-12-31 21:49:51
	% DurationCPUTime: 0.05s
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
	% StartTime: 2019-12-31 21:49:51
	% EndTime: 2019-12-31 21:49:51
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (68->10), mult. (36->12), div. (0->0), fcn. (18->6), ass. (0->13)
	t48 = qJD(1) + qJD(2);
	t55 = pkin(2) * t48;
	t54 = pkin(1) * qJD(1);
	t49 = qJ(1) + qJ(2);
	t47 = qJ(3) + t49;
	t42 = sin(t47);
	t43 = cos(t47);
	t44 = qJD(3) + t48;
	t53 = (-r_i_i_C(1) * t43 + r_i_i_C(2) * t42) * t44;
	t52 = (-r_i_i_C(1) * t42 - r_i_i_C(2) * t43) * t44;
	t51 = -cos(t49) * t55 + t53;
	t50 = -sin(t49) * t55 + t52;
	t1 = [-cos(qJ(1)) * t54 + t51, t51, t53, 0, 0; -sin(qJ(1)) * t54 + t50, t50, t52, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:49:51
	% EndTime: 2019-12-31 21:49:51
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (225->23), mult. (148->36), div. (0->0), fcn. (90->8), ass. (0->24)
	t73 = r_i_i_C(3) + pkin(8);
	t51 = qJ(1) + qJ(2);
	t49 = qJ(3) + t51;
	t44 = sin(t49);
	t45 = cos(t49);
	t53 = cos(qJ(4));
	t64 = qJD(4) * t53;
	t50 = qJD(1) + qJD(2);
	t46 = qJD(3) + t50;
	t52 = sin(qJ(4));
	t68 = t46 * t52;
	t72 = t44 * t68 - t45 * t64;
	t71 = t44 * t64 + t45 * t68;
	t70 = pkin(2) * t50;
	t67 = t46 * t53;
	t66 = pkin(1) * qJD(1);
	t65 = qJD(4) * t52;
	t61 = t44 * t65;
	t58 = t44 * t67 + t45 * t65;
	t57 = r_i_i_C(1) * t61 + ((-r_i_i_C(1) * t53 - pkin(3)) * t45 - t73 * t44) * t46 + t71 * r_i_i_C(2);
	t56 = -cos(t51) * t70 + t57;
	t55 = -t58 * r_i_i_C(1) + t72 * r_i_i_C(2) + (-pkin(3) * t44 + t45 * t73) * t46;
	t54 = -sin(t51) * t70 + t55;
	t1 = [-cos(qJ(1)) * t66 + t56, t56, t57, t72 * r_i_i_C(1) + t58 * r_i_i_C(2), 0; -sin(qJ(1)) * t66 + t54, t54, t55, (-t45 * t67 + t61) * r_i_i_C(2) - t71 * r_i_i_C(1), 0; 0, 0, 0, (-r_i_i_C(1) * t52 - r_i_i_C(2) * t53) * qJD(4), 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:49:52
	% EndTime: 2019-12-31 21:49:52
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (394->28), mult. (274->37), div. (0->0), fcn. (182->8), ass. (0->31)
	t181 = qJ(1) + qJ(2);
	t179 = qJ(3) + t181;
	t175 = cos(t179);
	t180 = qJD(1) + qJD(2);
	t176 = qJD(3) + t180;
	t200 = t175 * t176;
	t205 = pkin(4) + r_i_i_C(1);
	t206 = pkin(8) + r_i_i_C(2);
	t203 = r_i_i_C(3) + qJ(5);
	t183 = cos(qJ(4));
	t192 = t203 * t183;
	t182 = sin(qJ(4));
	t195 = t205 * t182;
	t190 = -t195 + t192;
	t204 = pkin(2) * t180;
	t202 = pkin(1) * qJD(1);
	t174 = sin(t179);
	t201 = t174 * t176;
	t199 = t176 * t182;
	t198 = qJD(4) * t182;
	t197 = qJD(4) * t183;
	t196 = t182 * qJD(5);
	t193 = t175 * t197;
	t191 = -t203 * t182 - t205 * t183;
	t189 = -pkin(3) + t191;
	t188 = t191 * qJD(4) + qJD(5) * t183;
	t187 = t189 * t201 + t206 * t200 + t203 * t193 + (-qJD(4) * t195 + t196) * t175;
	t186 = -sin(t181) * t204 + t187;
	t185 = t189 * t200 + (-qJD(4) * t192 - t206 * t176 + t205 * t198 - t196) * t174;
	t184 = -cos(t181) * t204 + t185;
	t1 = [-cos(qJ(1)) * t202 + t184, t184, t185, t188 * t175 - t190 * t201, -t174 * t199 + t193; -sin(qJ(1)) * t202 + t186, t186, t187, t188 * t174 + t190 * t200, t174 * t197 + t175 * t199; 0, 0, 0, t190 * qJD(4) + t196, t198;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end