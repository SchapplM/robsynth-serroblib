% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRRPP2
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 19:37
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RRRPP2_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP2_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP2_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRRPP2_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRPP2_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP2_jacobiaD_transl_sym_varpar: pkin has to be [7x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:37:02
	% EndTime: 2019-12-29 19:37:02
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:37:03
	% EndTime: 2019-12-29 19:37:03
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
	% StartTime: 2019-12-29 19:36:55
	% EndTime: 2019-12-29 19:36:55
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
	% StartTime: 2019-12-29 19:37:08
	% EndTime: 2019-12-29 19:37:08
	% DurationCPUTime: 0.18s
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
	t45 = -r_i_i_C(1) * t47 + t60 * r_i_i_C(2) + (-pkin(2) * t39 + t40 * t61) * t41;
	t1 = [-cos(qJ(1)) * t55 + t46, t46, t60 * r_i_i_C(1) + t47 * r_i_i_C(2), 0, 0; -sin(qJ(1)) * t55 + t45, t45, (-t40 * t56 + t50) * r_i_i_C(2) - t59 * r_i_i_C(1), 0, 0; 0, 0, (-r_i_i_C(1) * t43 - r_i_i_C(2) * t44) * qJD(3), 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:37:09
	% EndTime: 2019-12-29 19:37:10
	% DurationCPUTime: 0.33s
	% Computational Cost: add. (176->24), mult. (202->34), div. (0->0), fcn. (136->6), ass. (0->26)
	t173 = qJ(1) + qJ(2);
	t171 = cos(t173);
	t172 = qJD(1) + qJD(2);
	t190 = t171 * t172;
	t194 = pkin(3) + r_i_i_C(1);
	t195 = r_i_i_C(2) + pkin(7);
	t193 = r_i_i_C(3) + qJ(4);
	t175 = cos(qJ(3));
	t182 = t193 * t175;
	t174 = sin(qJ(3));
	t185 = t194 * t174;
	t180 = -t185 + t182;
	t192 = pkin(1) * qJD(1);
	t170 = sin(t173);
	t191 = t170 * t172;
	t189 = t172 * t174;
	t188 = qJD(3) * t174;
	t187 = qJD(3) * t175;
	t186 = t174 * qJD(4);
	t183 = t171 * t187;
	t181 = -t193 * t174 - t194 * t175;
	t179 = -pkin(2) + t181;
	t178 = t181 * qJD(3) + qJD(4) * t175;
	t177 = t179 * t191 + t195 * t190 + t193 * t183 + (-qJD(3) * t185 + t186) * t171;
	t176 = t179 * t190 + (-qJD(3) * t182 - t195 * t172 + t194 * t188 - t186) * t170;
	t1 = [-cos(qJ(1)) * t192 + t176, t176, t178 * t171 - t180 * t191, -t170 * t189 + t183, 0; -sin(qJ(1)) * t192 + t177, t177, t178 * t170 + t180 * t190, t170 * t187 + t171 * t189, 0; 0, 0, t180 * qJD(3) + t186, t188, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:37:08
	% EndTime: 2019-12-29 19:37:08
	% DurationCPUTime: 0.27s
	% Computational Cost: add. (233->29), mult. (254->37), div. (0->0), fcn. (171->6), ass. (0->26)
	t67 = r_i_i_C(2) + qJ(4);
	t73 = r_i_i_C(3) + qJ(5);
	t50 = sin(qJ(3));
	t62 = pkin(3) + pkin(4) + r_i_i_C(1);
	t58 = t62 * t50;
	t51 = cos(qJ(3));
	t59 = t67 * t51;
	t72 = -t59 + t58;
	t48 = qJD(1) + qJD(2);
	t57 = -t67 * t50 - t62 * t51;
	t71 = (-pkin(2) + t57) * t48 - qJD(5);
	t49 = qJ(1) + qJ(2);
	t46 = sin(t49);
	t70 = t48 * t46;
	t47 = cos(t49);
	t69 = t48 * t47;
	t68 = t48 * t50;
	t66 = pkin(1) * qJD(1);
	t65 = qJD(3) * t50;
	t64 = qJD(3) * t51;
	t63 = t50 * qJD(4);
	t60 = t47 * t64;
	t54 = t57 * qJD(3) + qJD(4) * t51;
	t53 = pkin(7) * t69 + t67 * t60 + t71 * t46 + (-qJD(3) * t58 - t73 * t48 + t63) * t47;
	t52 = t71 * t47 + t73 * t70 + (-pkin(7) * t48 - qJD(3) * t59 + t62 * t65 - t63) * t46;
	t1 = [-cos(qJ(1)) * t66 + t52, t52, t54 * t47 + t72 * t70, -t46 * t68 + t60, -t69; -sin(qJ(1)) * t66 + t53, t53, t54 * t46 - t69 * t72, t46 * t64 + t47 * t68, -t70; 0, 0, -qJD(3) * t72 + t63, t65, 0;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end