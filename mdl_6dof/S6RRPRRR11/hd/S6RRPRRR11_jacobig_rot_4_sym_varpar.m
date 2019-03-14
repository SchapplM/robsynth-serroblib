% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPRRR11
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:59
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRPRRR11_jacobig_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR11_jacobig_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR11_jacobig_rot_4_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:59:48
% EndTime: 2019-02-26 21:59:48
% DurationCPUTime: 0.01s
% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (7->4), ass. (0->4)
t59 = cos(qJ(1));
t58 = cos(qJ(2));
t57 = sin(qJ(1));
t1 = [0, t57, 0, t59 * t58, 0, 0; 0, -t59, 0, t57 * t58, 0, 0; 1, 0, 0, sin(qJ(2)) 0, 0;];
Jg_rot  = t1;