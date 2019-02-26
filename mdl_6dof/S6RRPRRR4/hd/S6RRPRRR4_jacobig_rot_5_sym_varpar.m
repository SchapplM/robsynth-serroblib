% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRRR4
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:55
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRPRRR4_jacobig_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR4_jacobig_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR4_jacobig_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:55:46
% EndTime: 2019-02-26 21:55:46
% DurationCPUTime: 0.02s
% Computational Cost: add. (15->5), mult. (42->12), div. (0->0), fcn. (65->8), ass. (0->15)
t109 = sin(pkin(12));
t111 = cos(pkin(12));
t113 = sin(qJ(2));
t115 = cos(qJ(2));
t117 = t109 * t113 - t111 * t115;
t116 = cos(qJ(1));
t114 = sin(qJ(1));
t112 = cos(pkin(6));
t110 = sin(pkin(6));
t108 = -t115 * t109 - t113 * t111;
t107 = t117 * t112;
t106 = t117 * t110;
t105 = -t114 * t107 - t116 * t108;
t104 = t116 * t107 - t114 * t108;
t1 = [0, t114 * t110, 0, t105, t105, 0; 0, -t116 * t110, 0, t104, t104, 0; 1, t112, 0, t106, t106, 0;];
Jg_rot  = t1;
