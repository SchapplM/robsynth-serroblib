% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRPP3
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:10
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRRRPP3_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP3_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP3_jacobig_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:10:01
% EndTime: 2019-02-26 20:10:01
% DurationCPUTime: 0.03s
% Computational Cost: add. (9->9), mult. (24->20), div. (0->0), fcn. (40->8), ass. (0->13)
t109 = sin(pkin(10));
t110 = sin(pkin(6));
t120 = t109 * t110;
t111 = cos(pkin(10));
t119 = t111 * t110;
t112 = cos(pkin(6));
t114 = sin(qJ(2));
t118 = t112 * t114;
t116 = cos(qJ(2));
t117 = t112 * t116;
t115 = cos(qJ(3));
t113 = sin(qJ(3));
t1 = [0, t120, t109 * t117 + t111 * t114 (-t109 * t118 + t111 * t116) * t113 - t115 * t120, 0, 0; 0, -t119, t109 * t114 - t111 * t117 (t109 * t116 + t111 * t118) * t113 + t115 * t119, 0, 0; 0, t112, -t110 * t116, t110 * t114 * t113 - t112 * t115, 0, 0;];
Jg_rot  = t1;
