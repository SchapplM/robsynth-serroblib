% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRRPP2
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:09
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRRRPP2_jacobigD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP2_jacobigD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPP2_jacobigD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP2_jacobigD_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:09:23
% EndTime: 2019-02-26 20:09:23
% DurationCPUTime: 0.07s
% Computational Cost: add. (12->10), mult. (48->29), div. (0->0), fcn. (50->8), ass. (0->16)
t153 = sin(pkin(6));
t156 = sin(qJ(3));
t166 = t153 * t156;
t157 = sin(qJ(2));
t165 = t153 * t157;
t155 = cos(pkin(6));
t164 = t155 * t157;
t159 = cos(qJ(2));
t163 = t155 * t159;
t162 = qJD(2) * t156;
t152 = sin(pkin(10));
t154 = cos(pkin(10));
t161 = t152 * t159 + t154 * t164;
t160 = -t152 * t164 + t154 * t159;
t158 = cos(qJ(3));
t1 = [0, 0, t160 * qJD(2) (t152 * t166 + t158 * t160) * qJD(3) + (-t152 * t163 - t154 * t157) * t162, 0, 0; 0, 0, t161 * qJD(2) (-t154 * t166 + t158 * t161) * qJD(3) + (-t152 * t157 + t154 * t163) * t162, 0, 0; 0, 0, qJD(2) * t165, t153 * t159 * t162 + (t155 * t156 + t158 * t165) * qJD(3), 0, 0;];
JgD_rot  = t1;
