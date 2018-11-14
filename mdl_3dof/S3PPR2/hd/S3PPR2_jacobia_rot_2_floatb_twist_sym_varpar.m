% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 2 (0=Basis) von
% S3PPR2
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
%
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
%
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d3]';
%
% Output:
% Ja_rot [3x3]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 10:11
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Ja_rot = S3PPR2_jacobia_rot_2_floatb_twist_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PPR2_jacobia_rot_2_floatb_twist_sym_varpar: qJ has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S3PPR2_jacobia_rot_2_floatb_twist_sym_varpar: pkin has to be [3x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_2_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:11:22
% EndTime: 2018-11-14 10:11:22
% DurationCPUTime: 0.01s
% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
t1 = [NaN, NaN, NaN; NaN, NaN, NaN; NaN, NaN, NaN;];
Ja_rot  = t1;
