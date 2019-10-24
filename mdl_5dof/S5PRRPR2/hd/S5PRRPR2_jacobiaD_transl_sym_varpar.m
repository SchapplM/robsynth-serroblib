% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5PRRPR2
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
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:30
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5PRRPR2_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR2_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR2_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5PRRPR2_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRPR2_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR2_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:30:16
	% EndTime: 2019-10-24 10:30:16
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:30:16
	% EndTime: 2019-10-24 10:30:16
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:30:16
	% EndTime: 2019-10-24 10:30:16
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (6->3), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->4)
	t32 = pkin(8) + qJ(2);
	t31 = cos(t32);
	t30 = sin(t32);
	t1 = [0, (-r_i_i_C(1) * t31 + r_i_i_C(2) * t30) * qJD(2), 0, 0, 0; 0, (-r_i_i_C(1) * t30 - r_i_i_C(2) * t31) * qJD(2), 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:30:16
	% EndTime: 2019-10-24 10:30:16
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (32->7), mult. (20->9), div. (0->0), fcn. (10->4), ass. (0->9)
	t45 = pkin(2) * qJD(2);
	t41 = pkin(8) + qJ(2);
	t40 = qJ(3) + t41;
	t38 = sin(t40);
	t39 = cos(t40);
	t42 = qJD(2) + qJD(3);
	t44 = (-r_i_i_C(1) * t39 + r_i_i_C(2) * t38) * t42;
	t43 = (-r_i_i_C(1) * t38 - r_i_i_C(2) * t39) * t42;
	t1 = [0, -cos(t41) * t45 + t44, t44, 0, 0; 0, -sin(t41) * t45 + t43, t43, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:30:16
	% EndTime: 2019-10-24 10:30:16
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (98->12), mult. (58->14), div. (0->0), fcn. (36->6), ass. (0->14)
	t38 = pkin(8) + qJ(2);
	t37 = qJ(3) + t38;
	t35 = sin(t37);
	t39 = qJD(2) + qJD(3);
	t48 = t39 * t35;
	t50 = qJ(4) + r_i_i_C(3);
	t49 = qJD(4) + r_i_i_C(2) * t39 * sin(pkin(9));
	t36 = cos(t37);
	t47 = t39 * t36;
	t46 = pkin(2) * qJD(2);
	t44 = -r_i_i_C(1) * cos(pkin(9)) - pkin(3);
	t43 = t49 * t35 + t44 * t48 + t50 * t47;
	t42 = -t50 * t48 + (t44 * t39 + t49) * t36;
	t1 = [0, -cos(t38) * t46 + t42, t42, t47, 0; 0, -sin(t38) * t46 + t43, t43, t48, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:30:16
	% EndTime: 2019-10-24 10:30:16
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (265->27), mult. (200->44), div. (0->0), fcn. (164->8), ass. (0->26)
	t191 = pkin(2) * qJD(2);
	t174 = pkin(8) + qJ(2);
	t173 = qJ(3) + t174;
	t171 = sin(t173);
	t175 = qJD(2) + qJD(3);
	t190 = t175 * t171;
	t172 = cos(t173);
	t189 = t175 * t172;
	t177 = cos(pkin(9));
	t178 = sin(qJ(5));
	t188 = t177 * t178;
	t179 = cos(qJ(5));
	t187 = t177 * t179;
	t186 = -t171 * t178 - t172 * t187;
	t185 = -t171 * t179 + t172 * t188;
	t184 = t171 * t187 - t172 * t178;
	t183 = t171 * t188 + t172 * t179;
	t176 = sin(pkin(9));
	t182 = -pkin(4) * t177 - pkin(3) + (-pkin(7) - r_i_i_C(3)) * t176;
	t164 = t186 * qJD(5) + t183 * t175;
	t165 = t185 * qJD(5) + t184 * t175;
	t181 = -t165 * r_i_i_C(1) + t164 * r_i_i_C(2) + qJ(4) * t189 + t171 * qJD(4) + t182 * t190;
	t166 = t184 * qJD(5) + t185 * t175;
	t167 = t183 * qJD(5) + t186 * t175;
	t180 = t166 * r_i_i_C(2) + t167 * r_i_i_C(1) + t172 * qJD(4) + (-qJ(4) * t171 + t182 * t172) * t175;
	t1 = [0, -cos(t174) * t191 + t180, t180, t189, t164 * r_i_i_C(1) + t165 * r_i_i_C(2); 0, -sin(t174) * t191 + t181, t181, t190, -t166 * r_i_i_C(1) + t167 * r_i_i_C(2); 0, 0, 0, 0, (-r_i_i_C(1) * t179 + r_i_i_C(2) * t178) * t176 * qJD(5);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end