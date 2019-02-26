% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRPPRR3
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta4]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:46
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRPPRR3_jacobigD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR3_jacobigD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR3_jacobigD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR3_jacobigD_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:45:50
% EndTime: 2019-02-26 19:45:50
% DurationCPUTime: 0.04s
% Computational Cost: add. (7->7), mult. (30->20), div. (0->0), fcn. (30->8), ass. (0->10)
t114 = cos(pkin(6));
t115 = sin(qJ(2));
t118 = t114 * t115;
t116 = cos(qJ(2));
t117 = t114 * t116;
t113 = cos(pkin(10));
t112 = cos(pkin(11));
t110 = sin(pkin(10));
t109 = sin(pkin(11));
t1 = [0, 0, 0, 0 ((-t110 * t117 - t113 * t115) * t109 - (-t110 * t118 + t113 * t116) * t112) * qJD(2), 0; 0, 0, 0, 0 ((-t110 * t115 + t113 * t117) * t109 - (t110 * t116 + t113 * t118) * t112) * qJD(2), 0; 0, 0, 0, 0 (t109 * t116 - t112 * t115) * sin(pkin(6)) * qJD(2), 0;];
JgD_rot  = t1;
