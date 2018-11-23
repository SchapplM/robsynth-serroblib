% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S5RRRRR1
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
%
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-16 14:52
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function JR_rot = S5RRRRR1_jacobiR_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR1_jacobiR_rot_4_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RRRRR1_jacobiR_rot_4_sym_varpar: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-16 14:52:29
% EndTime: 2018-11-16 14:52:29
% DurationCPUTime: 0.02s
% Computational Cost: add. (63->20), mult. (16->4), div. (0->0), fcn. (40->4), ass. (0->10)
t28 = qJ(2) + qJ(3) + qJ(4);
t26 = sin(t28);
t29 = sin(qJ(1));
t34 = t29 * t26;
t27 = cos(t28);
t33 = t29 * t27;
t30 = cos(qJ(1));
t32 = t30 * t26;
t31 = t30 * t27;
t1 = [-t33, -t32, -t32, -t32, 0; t31, -t34, -t34, -t34, 0; 0, -t27, -t27, -t27, 0; t34, -t31, -t31, -t31, 0; -t32, -t33, -t33, -t33, 0; 0, t26, t26, t26, 0; -t30, 0, 0, 0, 0; -t29, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
JR_rot  = t1;
