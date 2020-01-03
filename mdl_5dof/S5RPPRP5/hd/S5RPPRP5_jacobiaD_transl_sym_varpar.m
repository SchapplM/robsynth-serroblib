% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RPPRP5
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
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RPPRP5_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP5_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP5_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RPPRP5_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPRP5_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP5_jacobiaD_transl_sym_varpar: pkin has to be [7x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:54:09
	% EndTime: 2019-12-31 17:54:09
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:54:09
	% EndTime: 2019-12-31 17:54:09
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
	% StartTime: 2019-12-31 17:54:09
	% EndTime: 2019-12-31 17:54:09
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
	% StartTime: 2019-12-31 17:54:10
	% EndTime: 2019-12-31 17:54:10
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (17->12), mult. (48->17), div. (0->0), fcn. (34->4), ass. (0->9)
	t98 = r_i_i_C(2) + qJ(2);
	t92 = sin(qJ(1));
	t97 = qJD(1) * t92;
	t93 = cos(qJ(1));
	t96 = qJD(1) * t93;
	t90 = sin(pkin(7));
	t95 = qJD(3) * t90;
	t94 = -pkin(1) + (-pkin(2) - r_i_i_C(1)) * cos(pkin(7)) + (-r_i_i_C(3) - qJ(3)) * t90;
	t1 = [-t92 * t95 + t93 * qJD(2) + (-t98 * t92 + t94 * t93) * qJD(1), t96, -t90 * t97, 0, 0; t93 * t95 + t92 * qJD(2) + (t94 * t92 + t98 * t93) * qJD(1), t97, t90 * t96, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:54:09
	% EndTime: 2019-12-31 17:54:09
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (50->24), mult. (160->44), div. (0->0), fcn. (136->6), ass. (0->22)
	t36 = sin(qJ(1));
	t49 = qJD(1) * t36;
	t38 = cos(qJ(1));
	t48 = qJD(1) * t38;
	t33 = sin(pkin(7));
	t47 = qJD(3) * t33;
	t46 = pkin(6) + r_i_i_C(3) - qJ(2);
	t34 = cos(pkin(7));
	t35 = sin(qJ(4));
	t37 = cos(qJ(4));
	t45 = -t33 * t37 + t34 * t35;
	t44 = t33 * t35 + t34 * t37;
	t43 = t44 * t38;
	t42 = t44 * t36;
	t41 = qJD(1) * t45;
	t40 = t45 * qJD(4);
	t39 = -qJ(3) * t33 - pkin(1) + (-pkin(2) - pkin(3)) * t34;
	t32 = -qJD(1) * t43 + t36 * t40;
	t31 = qJD(4) * t42 + t38 * t41;
	t30 = qJD(1) * t42 + t38 * t40;
	t29 = -qJD(4) * t43 + t36 * t41;
	t1 = [-t36 * t47 + t32 * r_i_i_C(1) + t31 * r_i_i_C(2) + qJD(2) * t38 + (t46 * t36 + t39 * t38) * qJD(1), t48, -t33 * t49, r_i_i_C(1) * t29 + r_i_i_C(2) * t30, 0; t38 * t47 - t30 * r_i_i_C(1) + t29 * r_i_i_C(2) + t36 * qJD(2) + (t39 * t36 - t46 * t38) * qJD(1), t49, t33 * t48, -r_i_i_C(1) * t31 + r_i_i_C(2) * t32, 0; 0, 0, 0, (t45 * r_i_i_C(1) + t44 * r_i_i_C(2)) * qJD(4), 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:54:10
	% EndTime: 2019-12-31 17:54:10
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (103->34), mult. (325->49), div. (0->0), fcn. (296->6), ass. (0->29)
	t186 = sin(pkin(7));
	t190 = cos(qJ(4));
	t187 = cos(pkin(7));
	t188 = sin(qJ(4));
	t210 = t187 * t188;
	t217 = -t186 * t190 + t210;
	t197 = t186 * t188 + t187 * t190;
	t193 = t197 * qJD(4);
	t195 = t217 * qJD(5);
	t216 = -t186 * qJD(3) - t195;
	t215 = t217 * qJD(4);
	t214 = pkin(4) + r_i_i_C(1);
	t213 = r_i_i_C(3) + qJ(5);
	t189 = sin(qJ(1));
	t208 = qJD(1) * t189;
	t191 = cos(qJ(1));
	t207 = qJD(1) * t191;
	t204 = pkin(6) + r_i_i_C(2) - qJ(2);
	t203 = qJD(1) * t210;
	t202 = t186 * t208;
	t201 = t186 * t207;
	t196 = qJD(1) * t197;
	t194 = t197 * qJD(5);
	t192 = -qJ(3) * t186 - pkin(1) + (-pkin(2) - pkin(3)) * t187;
	t176 = -t215 * t189 + t191 * t196;
	t175 = t189 * t193 - t190 * t201 + t191 * t203;
	t174 = t189 * t196 + t215 * t191;
	t173 = t189 * t203 - t190 * t202 - t191 * t193;
	t1 = [t191 * qJD(2) - t214 * t176 - t213 * t175 + t216 * t189 + (t204 * t189 + t192 * t191) * qJD(1), t207, -t202, t214 * t173 - t213 * t174 + t191 * t194, -t173; t189 * qJD(2) - t214 * t174 - t213 * t173 - t216 * t191 + (t192 * t189 - t204 * t191) * qJD(1), t208, t201, -t214 * t175 + t213 * t176 + t189 * t194, t175; 0, 0, 0, -t213 * t193 + t214 * t215 - t195, -t215;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end