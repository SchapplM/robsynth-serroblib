% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S6PRPPRR3
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:46
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRPPRR3_jacobiaD_transl_3_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR3_jacobiaD_transl_3_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR3_jacobiaD_transl_3_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPPRR3_jacobiaD_transl_3_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR3_jacobiaD_transl_3_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:46:00
% EndTime: 2019-02-26 19:46:00
% DurationCPUTime: 0.04s
% Computational Cost: add. (24->13), mult. (82->26), div. (0->0), fcn. (72->6), ass. (0->16)
t132 = sin(pkin(10));
t134 = cos(pkin(10));
t137 = cos(qJ(2));
t135 = cos(pkin(6));
t136 = sin(qJ(2));
t142 = t135 * t136;
t146 = t132 * t142 - t134 * t137;
t145 = -pkin(2) - r_i_i_C(1);
t144 = r_i_i_C(3) + qJ(3);
t141 = t135 * t137;
t139 = qJD(2) * t144;
t138 = t132 * t137 + t134 * t142;
t133 = sin(pkin(6));
t130 = t146 * qJD(2);
t128 = t138 * qJD(2);
t1 = [0, -t146 * qJD(3) - t145 * t130 - (t132 * t141 + t134 * t136) * t139, -t130, 0, 0, 0; 0, t138 * qJD(3) + t145 * t128 - (t132 * t136 - t134 * t141) * t139, t128, 0, 0, 0; 0 (qJD(3) * t136 + (t145 * t136 + t144 * t137) * qJD(2)) * t133, t133 * qJD(2) * t136, 0, 0, 0;];
JaD_transl  = t1;
