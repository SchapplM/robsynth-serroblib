% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRRPP3
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
% Datum: 2019-12-29 19:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RRRPP3_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP3_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP3_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRRPP3_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRPP3_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP3_jacobiaD_transl_sym_varpar: pkin has to be [7x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:40:01
	% EndTime: 2019-12-29 19:40:01
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:40:01
	% EndTime: 2019-12-29 19:40:01
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
	% StartTime: 2019-12-29 19:40:01
	% EndTime: 2019-12-29 19:40:01
	% DurationCPUTime: 0.09s
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
	% StartTime: 2019-12-29 19:40:06
	% EndTime: 2019-12-29 19:40:07
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
	% StartTime: 2019-12-29 19:40:08
	% EndTime: 2019-12-29 19:40:08
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (176->25), mult. (202->36), div. (0->0), fcn. (136->6), ass. (0->24)
	t177 = sin(qJ(3));
	t178 = cos(qJ(3));
	t195 = r_i_i_C(3) + qJ(4);
	t196 = pkin(3) - r_i_i_C(2);
	t200 = -t195 * t177 - t196 * t178;
	t199 = pkin(7) + r_i_i_C(1);
	t175 = qJD(1) + qJD(2);
	t198 = (-pkin(2) + t200) * t175;
	t184 = t195 * t178;
	t183 = -t196 * t177 + t184;
	t194 = pkin(1) * qJD(1);
	t176 = qJ(1) + qJ(2);
	t174 = cos(t176);
	t193 = t174 * t175;
	t192 = t175 * t177;
	t191 = qJD(3) * t177;
	t190 = qJD(3) * t178;
	t189 = t177 * qJD(4);
	t186 = t174 * t190;
	t181 = t200 * qJD(3) + qJD(4) * t178;
	t173 = sin(t176);
	t180 = t198 * t173 + t195 * t186 + t199 * t193 + (-t191 * t196 + t189) * t174;
	t179 = t198 * t174 + (pkin(3) * t191 - t189 - t199 * t175 + (-r_i_i_C(2) * t177 - t184) * qJD(3)) * t173;
	t1 = [-cos(qJ(1)) * t194 + t179, t179, -t183 * t175 * t173 + t181 * t174, -t173 * t192 + t186, 0; -sin(qJ(1)) * t194 + t180, t180, t181 * t173 + t183 * t193, t173 * t190 + t174 * t192, 0; 0, 0, t183 * qJD(3) + t189, t191, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:40:02
	% EndTime: 2019-12-29 19:40:03
	% DurationCPUTime: 0.25s
	% Computational Cost: add. (240->29), mult. (270->43), div. (0->0), fcn. (187->6), ass. (0->28)
	t201 = r_i_i_C(2) + qJ(4);
	t181 = cos(qJ(3));
	t180 = sin(qJ(3));
	t193 = t180 * qJD(4);
	t205 = t181 * qJD(5) + t193;
	t191 = pkin(3) + r_i_i_C(3) + qJ(5);
	t188 = t191 * t180;
	t204 = -t201 * t181 + t188;
	t203 = pkin(7) + pkin(4) + r_i_i_C(1);
	t200 = pkin(1) * qJD(1);
	t179 = qJ(1) + qJ(2);
	t176 = sin(t179);
	t178 = qJD(1) + qJD(2);
	t199 = t176 * t178;
	t177 = cos(t179);
	t198 = t177 * t178;
	t197 = t178 * t180;
	t196 = t178 * t181;
	t195 = qJD(3) * t180;
	t194 = qJD(3) * t181;
	t190 = t176 * t195;
	t189 = t177 * t194;
	t187 = -t201 * t180 - t191 * t181;
	t185 = -pkin(2) + t187;
	t184 = t187 * qJD(3) + qJD(4) * t181 - qJD(5) * t180;
	t183 = t185 * t199 + t201 * t189 + t203 * t198 + (-qJD(3) * t188 + t205) * t177;
	t182 = (-t193 + (-t201 * qJD(3) - qJD(5)) * t181) * t176 + (-t203 * t176 + t185 * t177) * t178 + t191 * t190;
	t1 = [-cos(qJ(1)) * t200 + t182, t182, t184 * t177 + t204 * t199, -t176 * t197 + t189, -t176 * t196 - t177 * t195; -sin(qJ(1)) * t200 + t183, t183, t184 * t176 - t198 * t204, t176 * t194 + t177 * t197, t177 * t196 - t190; 0, 0, -qJD(3) * t204 + t205, t195, t194;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end