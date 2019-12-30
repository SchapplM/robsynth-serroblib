% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRRPP1
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 19:34
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RRRPP1_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP1_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP1_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRRPP1_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRPP1_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP1_jacobiaD_transl_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:33:47
	% EndTime: 2019-12-29 19:33:47
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:33:37
	% EndTime: 2019-12-29 19:33:37
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:33:42
	% EndTime: 2019-12-29 19:33:42
	% DurationCPUTime: 0.10s
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
	% StartTime: 2019-12-29 19:33:32
	% EndTime: 2019-12-29 19:33:32
	% DurationCPUTime: 0.19s
	% Computational Cost: add. (93->19), mult. (104->33), div. (0->0), fcn. (64->6), ass. (0->19)
	t61 = r_i_i_C(3) + pkin(7);
	t42 = qJ(1) + qJ(2);
	t39 = sin(t42);
	t40 = cos(t42);
	t44 = cos(qJ(3));
	t53 = qJD(3) * t44;
	t41 = qJD(1) + qJD(2);
	t43 = sin(qJ(3));
	t57 = t41 * t43;
	t60 = t39 * t57 - t40 * t53;
	t59 = t39 * t53 + t40 * t57;
	t56 = t41 * t44;
	t55 = pkin(1) * qJD(1);
	t54 = qJD(3) * t43;
	t50 = t39 * t54;
	t47 = t39 * t56 + t40 * t54;
	t46 = r_i_i_C(1) * t50 + ((-r_i_i_C(1) * t44 - pkin(2)) * t40 - t61 * t39) * t41 + t59 * r_i_i_C(2);
	t45 = -t47 * r_i_i_C(1) + t60 * r_i_i_C(2) + (-pkin(2) * t39 + t40 * t61) * t41;
	t1 = [-cos(qJ(1)) * t55 + t46, t46, t60 * r_i_i_C(1) + t47 * r_i_i_C(2), 0, 0; -sin(qJ(1)) * t55 + t45, t45, (-t40 * t56 + t50) * r_i_i_C(2) - t59 * r_i_i_C(1), 0, 0; 0, 0, (-r_i_i_C(1) * t43 - r_i_i_C(2) * t44) * qJD(3), 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:33:42
	% EndTime: 2019-12-29 19:33:42
	% DurationCPUTime: 0.28s
	% Computational Cost: add. (158->22), mult. (140->27), div. (0->0), fcn. (91->8), ass. (0->21)
	t52 = qJD(1) + qJD(2);
	t53 = qJ(3) + pkin(8);
	t48 = sin(t53);
	t49 = cos(t53);
	t62 = r_i_i_C(1) * t48 + r_i_i_C(2) * t49 + sin(qJ(3)) * pkin(3);
	t60 = t62 * qJD(3);
	t78 = -t60 - t52 * (-qJ(4) - pkin(7));
	t71 = r_i_i_C(2) * t48;
	t77 = t52 * t71 + qJD(4);
	t76 = -r_i_i_C(1) * t49 - cos(qJ(3)) * pkin(3);
	t54 = qJ(1) + qJ(2);
	t50 = sin(t54);
	t69 = t52 * t50;
	t51 = cos(t54);
	t68 = t52 * t51;
	t66 = pkin(1) * qJD(1);
	t63 = -pkin(2) + t76;
	t61 = qJD(3) * (t71 + t76);
	t59 = (t63 * t52 + t77) * t51 + (-r_i_i_C(3) * t52 - t78) * t50;
	t58 = r_i_i_C(3) * t68 + t77 * t50 + t78 * t51 + t63 * t69;
	t1 = [-cos(qJ(1)) * t66 + t59, t59, t51 * t61 + t62 * t69, t68, 0; -sin(qJ(1)) * t66 + t58, t58, t50 * t61 - t62 * t68, t69, 0; 0, 0, -t60, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:33:43
	% EndTime: 2019-12-29 19:33:44
	% DurationCPUTime: 0.29s
	% Computational Cost: add. (279->31), mult. (238->42), div. (0->0), fcn. (163->8), ass. (0->26)
	t182 = qJ(3) + pkin(8);
	t177 = sin(t182);
	t206 = pkin(4) + r_i_i_C(1);
	t208 = t206 * t177 + sin(qJ(3)) * pkin(3);
	t203 = r_i_i_C(3) + qJ(5);
	t178 = cos(t182);
	t194 = t203 * t178;
	t190 = -t208 + t194;
	t207 = -t203 * t177 - t206 * t178 - cos(qJ(3)) * pkin(3);
	t202 = pkin(1) * qJD(1);
	t183 = qJ(1) + qJ(2);
	t179 = sin(t183);
	t181 = qJD(1) + qJD(2);
	t201 = t181 * t179;
	t180 = cos(t183);
	t200 = t181 * t180;
	t199 = qJD(3) * t179;
	t198 = qJD(3) * t180;
	t197 = t177 * qJD(5);
	t195 = t178 * t198;
	t191 = -pkin(2) + t207;
	t189 = t207 * qJD(3) + qJD(5) * t178;
	t184 = -qJ(4) - pkin(7);
	t188 = t180 * t197 + r_i_i_C(2) * t200 + t179 * qJD(4) - t208 * t198 + (t191 * t179 - t180 * t184) * t181 + t203 * t195;
	t187 = t184 * t201 + t180 * qJD(4) + (-qJD(3) * t194 - t197) * t179 + (-r_i_i_C(2) * t179 + t191 * t180) * t181 + t208 * t199;
	t1 = [-cos(qJ(1)) * t202 + t187, t187, t189 * t180 - t190 * t201, t200, -t177 * t201 + t195; -sin(qJ(1)) * t202 + t188, t188, t189 * t179 + t190 * t200, t201, t177 * t200 + t178 * t199; 0, 0, t190 * qJD(3) + t197, 0, qJD(3) * t177;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end