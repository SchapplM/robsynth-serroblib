% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S4PPRP1
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,theta1]';
%
% Output:
% JaD_transl [3x4]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:19
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S4PPRP1_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRP1_jacobiaD_transl_4_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRP1_jacobiaD_transl_4_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S4PPRP1_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPRP1_jacobiaD_transl_4_sym_varpar: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:19:18
% EndTime: 2019-02-26 19:19:18
% DurationCPUTime: 0.03s
% Computational Cost: add. (20->9), mult. (54->12), div. (0->0), fcn. (48->4), ass. (0->11)
t21 = sin(pkin(5));
t22 = cos(pkin(5));
t23 = sin(qJ(3));
t24 = cos(qJ(3));
t30 = -t21 * t24 + t22 * t23;
t29 = pkin(3) + r_i_i_C(1);
t26 = -r_i_i_C(3) - qJ(4);
t25 = t21 * t23 + t22 * t24;
t19 = t25 * qJD(3);
t18 = t30 * qJD(3);
t1 = [0, 0, t25 * qJD(4) + t26 * t18 - t29 * t19, t19; 0, 0, -t30 * qJD(4) + t29 * t18 + t26 * t19, -t18; 0, 0, 0, 0;];
JaD_transl  = t1;
